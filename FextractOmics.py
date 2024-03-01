#!/usr/bin/env python
# coding: utf-8

# In[1]:


def open_protein(pdb_file,table):
    with open(pdb_file, 'r') as file:
        atoms = []
        seqaux = []
        seq = ''
        res = {}

        for line in file:
            if line.startswith('ATOM'):
                res[line[22:26]] = line[17:20]
                atom = [
                    line[6:11].strip(),    # Atom serial number
                    line[12:16].strip(),   # Atom name
                    line[17:20].strip(),   # Residue name
                    line[21].strip(),      # Chain identifier
                    line[22:26].strip(),   # Residue sequence number
                    line[30:38].strip(),   # X coordinate
                    line[38:46].strip(),   # Y coordinate
                    line[46:54].strip(),   # Z coordinate
                    line[76:78].strip()    # Element symbol
                ]
                atoms.append(atom)

        for aa in sorted(res):
            if res[aa] in table:
                seq += table[res[aa]]
            else:
                print(f"Skipping unrecognized amino acid: {res[aa]}")

        return atoms, seq


# In[2]:


def open_protein_ligand(mol2_file):
    count = 0
    mol2_file = open(mol2_file,'r')
    protein_atom = list()
    protein_bond = list()
    protein_tripos = {}
    for l in mol2_file:
        if l.startswith('@<TRIPOS>ATOM'):
            count += 1
        if l.startswith('@<TRIPOS>BOND'):
            count += 1
        if l.startswith('@<TRIPOS>UNITY_ATOM_ATTR'):
            continue 
        if len(l) > 20 and count==1:
            l = ' '.join(l.split())
            protein_atom.append (l.strip().split(' '))
            l = l.strip().split(' ')
            protein_tripos[l[0]] = l[0:6]
        if count==2:
            l = ' '.join(l.split())
            protein_bond.append (l.strip().split(' '))
        
    return protein_atom,protein_bond,protein_tripos


# In[3]:


def PP_Distance(protein,protein_ligand):
    import math
    import os

    #hydrogen_B = 0
    lines_total = 0
    lines_total2 = 0
    global refined
    refined = 0
    HC_total2 = 0.0
    VDW_total = 0.0
    repulsive = 0
    london = 0.0
    surface_tension = 0.0
    hydrophobicity = 0.0
    HC_allowed_1 = 0.0
    HC_allowed_2 = 0.0
        
    pldistmp = open("PL_DIST.tmp", "w")
    ppdistmp = open("PP_DIST.tmp", "w")

    #-----------------------CALCULOS----------------------------
    n = 0 #conta posicao do caractere
    i = 0

    #----------------------------------------- DADOS DO LIGANTE ------------------------------
    atomo = ""
    nome_atomo = ""
    x = ""
    y = ""
    z = ""
    enter = ""
    num_atomo = 0
    num_atomo_copia = 0
    num_x = 0.0
    num_y = 0.0
    num_z = 0.0

    # Distancia permitida
    permitido = 0.0
    
    # RADII 1 E RADII 2
    r1 = 0.0
    r2 = 0.0

    #------------------------------------------ DADOS DA PROTEINA----------------------------------
    p_atomo = ""
    p_nome_atomo = ""
    p_aminoacido = ""
    p_num_aminoacido = ""
    p_x = ""
    p_y = ""
    p_z = ""
    chain = ""
    num_p_atomo = 0
    num_p_atomo_copia = 0
    num_px = 0.0
    num_py = 0.0
    num_pz = 0.0
    counter = 0
    
    RT_found = open("RT_found.tmp", "w")

    #------------------------------------------------------------------------------------------
    line = "*"
    hibr = ""
    v = []
    found_hibr = 0
    found_p = 0
    
    ldist = ''
    pdist = ''
    
    i=0
    for linel in protein_ligand:
        found_hibr = 0        
        v.clear()
       
        atomo = linel[0]
        num_atomo = linel[0]

        nome_atomo = linel[1]
        num_x = float(linel[2])
        x = linel[2]
        num_y = float(linel[3])
        y = linel[3]

        num_z = float(linel[4])
        z = linel[4]

        hibr = linel[5]
        l_num_aminoacido = linel[6]

        # Ver isso aqui com Raquel
        if hibr.startswith("N.3"):
            r2 = 1.87
            found_hibr = 1
        elif hibr.startswith("N.1") or hibr.startswith("N.2") or hibr.startswith("N.ar") or hibr.startswith("N.pl3"):
            r2 = 1.86
            found_hibr = 1
        elif hibr.startswith("N.am"):
            r2 = 1.83
            found_hibr = 1
        elif hibr.startswith("C.3"):
            r2 = 1.94
            found_hibr = 1
        elif hibr.startswith("C.1") or hibr.startswith("C.2"):
            r2 = 1.90
            found_hibr = 1
        elif hibr.startswith("C.ar"):
            r2 = 1.85
            found_hibr = 1
        elif hibr.startswith("O.3"):
            r2 = 1.74
            found_hibr = 1
        elif hibr.startswith("O.2") or hibr.startswith("O.co2"):
            r2 = 1.66
            found_hibr = 1
        elif hibr.startswith("F") or hibr.startswith("O.w"):
            r2 = 1.77
            found_hibr = 1
        elif hibr.startswith("S"):
            r2 = 2.01
            found_hibr = 1
        elif hibr.startswith("S.3"):
            r2 = 2.09
            found_hibr = 1
        elif hibr.startswith("Cl"):
            r2 = 2.00
            found_hibr = 1
        elif hibr.startswith("Br"):
            r2 = 2.22
            found_hibr = 1
        elif hibr.startswith("I"):
            r2 = 2.42
            found_hibr = 1
        elif hibr.startswith("P"):
            r2 = 2.03
            found_hibr = 1

        #--------LIGANTE-------------RADII-Raio-de-VDW-PERMITIDOS
        DIS2 = 0
        j=0
        for linep in protein:
            if(j<=i):
                j+=1
                continue
            found_p = 0

            if linep[:1] == "" or linep[:1] == "\n" or linep[:1] == "\0":
                break
            p_atomo = linep[0]
            num_p_atomo = int(linep[0])
            p_nome_atomo = linep[1]
            p_aminoacido = linep[2]
            chain = linep[3]
            p_num_aminoacido = linep[4]
            num_px = float(linep[5])
            p_x = linep[5]
            num_py = float(linep[6])
            p_y = linep[6]
            num_pz = float(linep[7])
            p_z = linep[7]
            
            if p_num_aminoacido == l_num_aminoacido:
                j+=1
                continue

            DX = (num_x - num_px)
            DY = (num_y - num_py)
            DZ = (num_z - num_pz)
            EDX = (DX * DX)
            EDY = (DY * DY)
            EDZ = (DZ * DZ)
            SD = (EDX + EDY + EDZ)
            
            DIS2 = math.sqrt(SD)
            
            if DIS2 == 0.0:
                DIS2 = 0.000001
            
            # DIS2 = DIS
            # CONTINUACAO ACHA DISTANCIA PERMITIDA
            if p_nome_atomo == "N":
                r1 = 1.65; found_p = 1
            if p_nome_atomo == "O":
                r1 = 1.40; found_p = 1
            if p_nome_atomo == "NZ":
                r1 = 1.50; found_p = 1
            if p_nome_atomo == "S":
                r1 = 1.85;found_p = 1;
            if p_nome_atomo.startswith('C'):
                r1 = 1.76
                found_p = 2

            if p_nome_atomo.startswith(('CA', 'CB', 'CD', 'CE', 'CG', 'CH', 'CZ')):
                r1 = 1.87
                found_p = 2

            '''
            1.65    N, NE, NH1, NH2, ND2, NE2, ND1
            1.40    O, OD1, OD2, OE1, OE2, OH, 
            1.85    SG
            1.50    NZ
            '''

            permitido = r1 + r2
            permitido1 = permitido - 0.7
            extended = 0.7

            # RT
            if DIS2 <= permitido + extended and not nome_atomo.startswith('H') and found_hibr == 1 and found_p == 1:
                if atomo != '':
                    RT_found.write(atomo + '\n')

            # VDW
            # [JC] Van der Waals interactions (VDWs)
            if not hibr.startswith('H') and not atomo.startswith('H') and found_hibr == 1 and found_p != 0:
                # function for VDW
                VDW_radii = r1 + r2

                var_0 = VDW_radii / DIS2
                exp1 = 8.0
                exp2 = 4.0
                var_1 = var_0 ** exp1
                var_2 = var_0 ** exp2
                VDW_int = var_1 - (2 * var_2)

                if VDW_int >= 100:
                    VDW_int = 100
                if VDW_int <= 100:
                    VDW_total += VDW_int
                VDW_int = 0

            # HC
            # [JC] Hydrophobic contacts (HCs)
            HC = 0
            HC2 = 0
            HC_exp = 0

            if p_nome_atomo.startswith('C') and nome_atomo.startswith('C') and found_hibr == 1 and found_p != 0:
                HC_VDW = permitido
                HC_allowed_1 = HC_VDW + 0.5
                HC_allowed_2 = HC_VDW + 2.0

            # ----------------------------------------- HC2 -----------------------------------------
            if DIS2 <= HC_allowed_1:
                HC2 = 1

            if DIS2 > HC_allowed_1 and DIS2 <= HC_allowed_2:
                HC2 = (1 / 1.5) * ((HC_VDW + 2.0) ** 2 - DIS2 ** 2)

            if DIS2 > HC_allowed_2:
                HC2 = 0

            HC_total2 += HC2
            if DIS2 <= permitido and nome_atomo[:1] != "H" and found_hibr == 1:
                repulsive += 1
            if london + 1 / DIS2 != float("inf"):
                london += 1 / DIS2

            # ----------------------------------------- HC2 ------------------------------------------
            # [JC] Hydrophobic contacts (HCs) 2
            pldistmp_aux = ""
            ppdistmp_aux = ""
            if DIS2 >= permitido1 and DIS2 <= permitido and found_hibr == 1 and found_p == 1:
                if (p_nome_atomo[:1] == "N" or p_nome_atomo[:1] == "O" or p_nome_atomo[:1] == "S") and (nome_atomo[:1] == "N" or nome_atomo[:1] == "O" or nome_atomo[:1] == "S"):
                    n = 0

                    pldistmp_aux += atomo + " " + nome_atomo
                    pldistmp_aux += " "

                    pldistmp_aux += x
                    pldistmp_aux += " "

                    pldistmp_aux += y
                    pldistmp_aux += " "

                    pldistmp_aux += z + " " + hibr + "\n"
                    pldistmp.write(pldistmp_aux)
                    
                    # proteina: salva arquivo tmp com os atomos mais proximos
                    n = 0

                    ppdistmp_aux += p_atomo + " " + p_nome_atomo
                    ppdistmp_aux += " "
                    ppdistmp_aux += p_aminoacido
                    ppdistmp_aux += " "

                    ppdistmp_aux += chain
                    ppdistmp_aux += " "

                    ppdistmp_aux += p_num_aminoacido
                    ppdistmp_aux += " "

                    ppdistmp_aux += p_x
                    ppdistmp_aux += " "

                    ppdistmp_aux += p_y
                    ppdistmp_aux += " "

                    ppdistmp_aux += p_z + "\n"
                    ppdistmp.write(ppdistmp_aux)
                    n = 0
            j+=1
        i+=1
    aux = [HC_total2, VDW_total, 0.0, repulsive, london]
    return ('PP_DIST.tmp',aux)


# In[4]:


def raiz_proteina(PP_DIST,protein):
    root = open("P_ROOT_0.tmp", "w")
    AD = open(PP_DIST, "r")
    p_atomo = ""
    p_nome_atomo = ""
    p_aminoacido = ""
    p_num_aminoacido = ""
    px = ""
    py = ""
    pz = ""
    p_enter = ""
    num_p_num_aminoacido = 0
    num_pp_num_aminoacido = 0
    num_p_atomo = 0
    num_p_atomo_copia = 0
    num_px = 0.0
    num_py = 0.0
    num_pz = 0.0
    n = 0
    pp_atomo = ""
    pp_nome_atomo = ""
    pp_aminoacido = ""
    pp_num_aminoacido = ""
    ppx = ""
    ppy = ""
    ppz = ""
    pp_enter = ""
    num_pp_atomo = 0
    num_pp_atomo_copia = 0
    num_ppx = 0.0
    num_ppy = 0.0
    num_ppz = 0.0
    # tentar abrir arquivo da proteina
    for l in AD:
        l = l.strip().split(' ')
        pp_atomo = l[0]
        pp_nome_atomo = l[1]
        pp_aminoacido = l[2]
        pp_num_aminoacido = l[4]
        ppx = float(l[5])
        ppy = float(l[6])
        ppz = float(l[7])
        pp_enter = l[3]
        #AD_root = open("protein.tmp")
        aa = 0
        
        for lp in protein:
            aa += 1
            p_atomo = lp[0] #AD_root.read(7)
            p_nome_atomo = lp[1] #AD_root.read(4)
            p_aminoacido = lp[2] #AD_root.read(3)
            p_num_aminoacido = lp[4] #AD_root.read(7)
            px = lp[5] #AD_root.read(8)
            py = lp[6] #AD_root.read(8)
            pz = lp[7] #AD_root.read(8)
            p_enter = lp[3] #AD_root.read(2)
                        
            # -------------------------------N-PROLINA
            # [JC]
            if pp_nome_atomo[0]=='N' and pp_aminoacido.startswith('PRO'):
                if aa!=1:
                    num_p_num_aminoacido = int(p_num_aminoacido)
                    num_pp_num_aminoacido = int(pp_num_aminoacido)
                    if p_nome_atomo=='C' and num_p_num_aminoacido==(num_pp_num_aminoacido -1 ):
                        # cria arquivo temp com as raizes
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                    if p_nome_atomo.startswith('CA') and pp_aminoacido==p_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                    if p_nome_atomo.startswith('CD') and pp_aminoacido == p_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
            
                if aa == 1:  # se não for o primeiro aminoácido
                    if p_nome_atomo.startswith('CA') and pp_aminoacido == p_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # N-TODAS-PROTEINAS-EXCETO-PRO
            if pp_nome_atomo == 'N' and not pp_aminoacido.startswith('PRO'):
                # N
                if aa != 1:  # se nao for o prmeiro aminoacido
                    num_p_num_aminoacido = int(p_num_aminoacido)
                    num_pp_num_aminoacido = int(pp_num_aminoacido)
                    num_pp_num_aminoacido_2 = num_pp_num_aminoacido - 1
                    if p_nome_atomo == 'C' and num_p_num_aminoacido == num_pp_num_aminoacido_2:
                        # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter << endl;
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                    # num_p_num_aminoacido == (num_pp_num_aminoacido - 1))
                    if p_nome_atomo.startswith('CA') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                # if aa!=1
                if aa == 1:  # se for o primeiro aminoáido
                    if p_nome_atomo.startswith('CA') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ----------------------------------NZ-LYS
            if pp_nome_atomo.startswith('NZ') and pp_aminoacido.startswith('LYS'):
                if p_nome_atomo.startswith('CE') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            #----------------------------------NE-ARG
            if pp_nome_atomo.startswith('NE') and pp_aminoacido.startswith('ARG'):
                if p_nome_atomo.startswith('CD') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                if p_nome_atomo.startswith('CZ') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ----------------------------------NH1/2-ARG
            if pp_nome_atomo.startswith('NH1') and pp_aminoacido.startswith('ARG'):
                if p_nome_atomo.startswith('CZ') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    # cria arquivo temp com as raizes
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            #----------------------------------ND1-HIS
            if pp_nome_atomo.startswith('ND1') and pp_aminoacido.startswith('HIS'):
                if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                if p_nome_atomo.startswith('CE1') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ----------------------------------NE2-HIS
            if pp_nome_atomo.startswith('NE2') and pp_aminoacido.startswith('HIS'):
                if p_nome_atomo.startswith('CD2') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                if p_nome_atomo.startswith('CE1') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ---------------------------ND2--ASN
            if pp_nome_atomo.startswith('ND2') and pp_aminoacido.startswith('ASN'):
                if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ---------------------------NE1--GLN
            # if (pp_nome_atomo.[1]=='E' and pp_nome_atomo[2]=='2' and pp_aminoacido[0]=='G' and pp_aminoacido[1]=='L' and pp_aminoacido[2]=='N'):
            if pp_nome_atomo.startswith('NE2') and pp_aminoacido.startswith('GLN'):
                if p_nome_atomo.startswith('CD') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            # ---------------------------NE1--TRP                   
            if pp_nome_atomo.startswith('NE1') and pp_aminoacido.startswith('TRP'):
                if p_nome_atomo.startswith('CD1') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                if p_nome_atomo.startswith('CE2') and pp_aminoacido.startswith(p_aminoacido[0:3] and pp_num_aminoacido == p_num_aminoacido):
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                     
            if pp_nome_atomo[0] == 'O':
                # ----------------------------------O
                if pp_nome_atomo == 'O' and p_nome_atomo == 'C' and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << px << "*" << py << "*" << pz << "*" << p_enter;
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                # ---------------------------OD1/2--ASP
                if (pp_nome_atomo.startswith('OD1') or pp_nome_atomo.startswith('OD2')) and pp_aminoacido.startswith('ASP'):
                    if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                if (pp_nome_atomo.startswith('OE1') or pp_nome_atomo.startswith('OE2')) and pp_aminoacido.startswith('GLU'):
                    if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                # ---------------------------OD1--ASN    
                if pp_nome_atomo.startswith('OD1') and pp_aminoacido.startswith('ASN'):
                    if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                # ---------------------------OG--SER    
                if pp_nome_atomo.startswith('OG') and pp_aminoacido.startswith('SER'):
                    if p_nome_atomo.startswith('CB') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                #---------------------------OG1--THR	
                if pp_nome_atomo.startswith('OG1') and pp_aminoacido.startswith('THR'):
                    if p_nome_atomo.startswith('CB') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                # ---------------------------OH--TYR
                if pp_nome_atomo.startswith('OH') and pp_aminoacido.startswith('TYR'):
                    if p_nome_atomo.startswith('CZ') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

                # ---------------------------OE1--GLN
                if pp_nome_atomo.startswith('OE1') and pp_aminoacido.startswith('GLN'):
                    if p_nome_atomo.startswith('CD') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                # ---------------------------OXT--TODOS		
                if pp_nome_atomo.startswith('OXT'):
                    if p_nome_atomo[0] == 'C' and  pp_aminoacido.startswith(p_aminoacido[0:2]) and pp_num_aminoacido == p_num_aminoacido:
                        root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
            # -------------------------------------S
            # --------------------------------SD-MET
            if pp_nome_atomo.startswith('SD') and pp_aminoacido.startswith('MET'):
                if p_nome_atomo.startswith('CE') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')
                if p_nome_atomo.startswith('CG') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

            if pp_nome_atomo.startswith('SG') and pp_aminoacido.startswith('CYS'):
                if p_nome_atomo.startswith('CB') and pp_aminoacido.startswith(p_aminoacido[0:3]) and pp_num_aminoacido == p_num_aminoacido:
                    root.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + px + ' ' + py + ' ' + pz + ' ' + p_enter +'\n')

        root.write("B"+' '+pp_atomo+' '+pp_nome_atomo+' '+pp_aminoacido+' '+pp_num_aminoacido+' '+str(ppx)+' '+str(ppy)+' '+str(ppz)+' '+pp_enter+'\n')
    AD.close();
    root.close();


# In[5]:


def P_centro_geom():
    
    # r1 = raiz1
    pos = 0.0
    conta_r = 0.0
    tipo = ''
    p_atomo = ''
    p_nome_atomo = ''
    p_aminoacido = ''
    p_num_aminoacido = ''
    px = 0.0
    py = 0.0
    pz = 0.0
    p_enter = ''
    
    # raiz 1 
    r1_atomo = ''
    r1_nome_atomo = ''
    r1_aminoacido = ''
    r1_num_aminoacido = ''
    r1x = 0.0
    r1y = 0.0
    r1z = 0.0
    r1_enter = ''
    # raiz 2 
    r2_atomo = ''
    r2_nome_atomo = ''
    r2_aminoacido = ''
    r2_num_aminoacido = ''
    r2x = 0.0
    r2y = 0.0
    r2z = 0.0
    r2_enter = ''
    
    #--------converte caracteres para 2 raizes
    num_r1x = 0.0
    num_r1y = 0.0
    num_r1z = 0.0
    num_r2x = 0.0
    num_r2y = 0.0
    num_r2z = 0.0
    
    #--------converte caracteres para 3 raizes
    num_r3x = 0.0
    num_r3y = 0.0
    num_r3z = 0.0
    
    # raiz 3 
    r3_atomo = 0.0
    r3_nome_atomo = ''
    r3_aminoacido = ''
    r3_num_aminoacido = 0.0
    r3x = 0.0
    r3y = 0.0
    r3z = 0.0
    r3_enter = ''

    #--------Valores para calcular o centro geometrico (geometric center)
    gcx2 = 0.0
    gcy2 = 0.0
    gcz2 = 0.0
    gcx2t = 0.0
    gcy2t = 0.0
    gcz2t = 0.0
    gcx3 = 0.0
    gcy3 = 0.0
    gcz3 = 0.0
    gcx3t = 0.0
    gcy3t = 0.0
    gcz3t = 0.0
    
    p_root_f = open('P_ROOT.tmp', 'w')  # final
       
    with open("P_ROOT_0.tmp", "r") as p_root:
        if not p_root:
            print("Could not create temp file!")

        for l in p_root:
            if len (l) > 5:
                l = l.strip().split(' ')
                #   0     1     2     3     4       5         6       7      8    
                #   R    308    CA   GLY   21     3.307     10.759   7.547   B
                # pega uma linha
                tipo = l[0]
                p_atomo = l[1]
                p_nome_atomo = l[2] 
                p_aminoacido = l[3] 
                p_num_aminoacido = l[4] 
                px = l[5]
                py = l[6]
                pz = l[7]
                p_enter = l[8]

                # se nao for raiz e tiver so uma raiz
                # ... o resto do código aqui ...
                if tipo == 'B' and conta_r == 1:
                    conta_r = 0
                    p_root_f.write("R" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + str(r1x) + ' ' + str(r1y) + ' ' + str(r1z) + ' ' + p_enter + '\n')
                    p_root_f.write("B" + ' ' + p_atomo + ' ' + p_nome_atomo + ' ' + p_aminoacido + ' ' + p_num_aminoacido + ' ' + str(px) + ' ' + str(py) + ' ' + str(pz) + ' ' + p_enter + '\n')
                    continue # goto get # Ver isso

                # se a primeira linha for uma raiz
                if tipo == 'R' and conta_r == 0:
                    conta_r += 1
                    r1x = px
                    r1y = py
                    r1z = pz
                    #atualizar r1x, r1y,r1z
                    continue

                # se a primeira linha for uma raiz  
                if tipo == 'R' and conta_r == 1:
                    conta_r += 1
                    #atualizar r2x, r2y,r2z
                    r2x = px
                    r2y = py
                    r2z = pz
                    continue

                if tipo == 'R' and conta_r == 2:
                    conta_r += 1
                    #atualizar r3x, r3y,r3z
                    r3x = px
                    r3y = py
                    r3z = pz
                    continue

                if tipo == 'B' and (conta_r == 2 or conta_r == 3):
                    num_r1x = float(r1x)
                    num_r1y = float(r1y)
                    num_r1z = float(r1z)

                    num_r2x = float(r2x)
                    num_r2y = float(r2y)
                    num_r2z = float(r2z)

                    if conta_r == 2:
                        gcx2t = num_r1x + num_r2x
                        gcx2 = gcx2t / 2
                        gcy2t = num_r1y + num_r2y
                        gcy2 = gcy2t / 2
                        gcz2t = num_r1z + num_r2z
                        gcz2 = gcz2t / 2

                        if -99.9999 <= gcx2 <= 999.9999:
                            letrasx = gcx2
                        if -999.9999 <= gcx2 <= -100.0000:
                            letrasx = gcx2
                        if -99.9999 <= gcy2 <= 999.9999:
                            letrasy = gcy2
                        if -999.9999 <= gcy2 <= -100.0000:
                            letrasy = gcy2
                        if -99.9999 <= gcz2 <= 999.9999:
                            letrasz = gcz2
                        if -999.9999 <= gcz2 <= -100.0000:
                            letrasz = gcz2

                        p_root_f.write("R geom_center(2) " +" "+ str(letrasx) +" "+ str(letrasy) +" "+ str(letrasz) +" "+ "\n")
                        conta_r = 0

                    if conta_r == 3:

                        num_r3x = float(r3x)
                        num_r3y = float(r3y)
                        num_r3z = float(r3z)
                        gcx3t = (num_r1x + num_r2x + num_r3x)
                        gcx3 = (gcx3t / 3)
                        gcy3t = (num_r1y + num_r2y + num_r3y)
                        gcy3 = (gcy3t / 3)
                        gcz3t = (num_r1z + num_r2z + num_r3z)
                        gcz3 = (gcz3t / 3)

                        if -99.9999 <= gcx3 <= 999.9999:
                            letrasx = gcx3
                        if -999.9999 <= gcx3 <= -100.0000:
                            letrasx = gcx3
                        if -99.9999 <= gcy3 <= 999.9999:
                            letrasy = gcy3
                        if -999.9999 <= gcy3 <= -100.0000:
                            letrasy = gcy3
                        if -99.9999 <= gcz3 <= 999.9999:
                            letrasz = gcz3
                        if -999.9999 <= gcz3 <= -100.0000:
                            letrasz = gcz3


                        p_root_f.write("R geom_center(3) " +" "+ str(letrasx) +" "+ str(letrasy) +" "+ str(letrasz) +" "+ "\n")
                        conta_r = 0
                p_root_f.write(tipo +" "+ p_atomo +" "+ p_nome_atomo +" "+ p_aminoacido +" "+ p_num_aminoacido +" "+ px +" "+ py +" "+ pz +" "+ p_enter + "\n")
        if tipo != 'R' and tipo != 'B' and tipo != '\0':
            print("ERRO! |" + tipo + "|")

    p_root.close() 
    p_root_f.close()


# In[6]:


def L_centro_geom():

    # r1 = raiz1
    pos = 0.0
    conta_r = 0.0
    tipo = ''
    l_atomo = ''
    l_nome_atomo = ''
    l_aminoacido = ''
    l_num_aminoacido = ''
    lx = 0.0
    ly = 0.0
    lz = 0.0
    l_enter = ''

    # raiz 1 
    r1_atomo = ''
    r1_nome_atomo = ''
    r1_aminoacido = ''
    r1_num_aminoacido = ''
    r1x = 0.0
    r1y = 0.0
    r1z = 0.0
    r1_enter = ''
    # raiz 2 
    r2_atomo = ''
    r2_nome_atomo = ''
    r2_aminoacido = ''
    r2_num_aminoacido = ''
    r2x = 0.0
    r2y = 0.0
    r2z = 0.0
    r2_enter = ''
    
    #--------converte caracteres para 2 raizes
    num_r1x = 0.0
    num_r1y = 0.0
    num_r1z = 0.0
    num_r2x = 0.0
    num_r2y = 0.0
    num_r2z = 0.0
    
    #--------converte caracteres para 3 raizes
    num_r3x = 0.0
    num_r3y = 0.0
    num_r3z = 0.0
    
    # raiz 3 
    r3_atomo = 0.0
    r3_nome_atomo = ''
    r3_aminoacido = ''
    r3_num_aminoacido = 0.0
    r3x = 0.0
    r3y = 0.0
    r3z = 0.0
    r3_enter = ''

    #--------Valores para calcular o centro geometrico (geometric center)
    gcx2 = 0.0
    gcy2 = 0.0
    gcz2 = 0.0
    gcx2t = 0.0
    gcy2t = 0.0
    gcz2t = 0.0
    gcx3 = 0.0
    gcy3 = 0.0
    gcz3 = 0.0
    gcx3t = 0.0
    gcy3t = 0.0
    gcz3t = 0.0
    

    #--------valores para calcular o centro geometrico (geometric center)
    gcx2 = 0.0
    gcy2 = 0.0
    gcz2 = 0.0
    gcx2t = 0.0
    gcy2t = 0.0
    gcz2t = 0.0
    gcx3 = 0.0
    gcy3 = 0.0
    gcz3 = 0.0
    gcx3t = 0.0
    gcy3t = 0.0
    gcz3t = 0.0
    
    l_root_f = open("L_ROOT.tmp", "w")
    l_root = open("L_ROOT_0.tmp", "r")
    
    with open("L_ROOT_0.tmp", "r") as p_root:
        if not l_root:
            print("Could not create temp file!")

        for l in p_root:
            if len (l) > 5:
                l = l.strip().split(' ')
                #   0     1     2     3     4       5         6       7      8    
                #   R    308    CA   GLY   21     3.307     10.759   7.547   B
                #        12      N                14.5670   1.7920   7.9420  N.am
                #   B    4       O                16.5200    8.5150   25.9400 O.co2
                # pega uma linha
                if(l[1].startswith('geom_center')):
                    tipo = l[0]
                    l_atomo = l[1]
                    lx = l[2]
                    ly = l[3]
                    lz = l[4]
                    l_nome_atomo = ''
                    l_enter = ''
                else:
                    tipo = l[0]
                    l_atomo = l[1]
                    l_nome_atomo = l[2] 
                    lx = l[3]
                    ly = l[4]
                    lz = l[5]
                    l_enter = l[6]

                # se nao for raiz e tiver so uma raiz
                # ... o resto do código aqui ...
                if tipo == 'B' and conta_r == 1:
                    conta_r = 0
                    l_root_f.write("R" + ' ' + l_atomo + ' ' + l_nome_atomo + ' ' + str(r1x) + ' ' + str(r1y) + ' ' + str(r1z) + ' ' + l_enter + '\n')
                    l_root_f.write("B" + ' ' + l_atomo + ' ' + l_nome_atomo + ' ' + str(lx) + ' '  + str(ly)  + ' ' + str(lz)  + ' ' + l_enter + '\n')
                    continue # goto get # Ver isso

                # se a primeira linha for uma raiz
                if tipo == 'R' and conta_r == 0:
                    conta_r += 1
                    r1x = lx
                    r1y = ly
                    r1z = lz
                    continue

                # se a primeira linha for uma raiz  
                if tipo == 'R' and conta_r == 1:
                    conta_r += 1
                    r2x = lx
                    r2y = ly
                    r2z = lz
                    continue

                if tipo == 'R' and conta_r == 2:
                    conta_r += 1
                    r3x = lx
                    r3y = ly
                    r3z = lz
                    continue

                if tipo == 'B' and (conta_r == 2 or conta_r == 3):
                    num_r1x = float(r1x)
                    num_r1y = float(r1y)
                    num_r1z = float(r1z)

                    num_r2x = float(r2x)
                    num_r2y = float(r2y)
                    num_r2z = float(r2z)

                    if conta_r == 2:
                        gcx2t = num_r1x + num_r2x
                        gcx2 = gcx2t / 2
                        gcy2t = num_r1y + num_r2y
                        gcy2 = gcy2t / 2
                        gcz2t = num_r1z + num_r2z
                        gcz2 = gcz2t / 2

                        if -99.9999 <= gcx2 <= 999.9999:
                            letrasx = gcx2
                        if -999.9999 <= gcx2 <= -100.0000:
                            letrasx = gcx2
                        if -99.9999 <= gcy2 <= 999.9999:
                            letrasy = gcy2
                        if -999.9999 <= gcy2 <= -100.0000:
                            letrasy = gcy2
                        if -99.9999 <= gcz2 <= 999.9999:
                            letrasz = gcz2
                        if -999.9999 <= gcz2 <= -100.0000:
                            letrasz = gcz

                        l_root_f.write("R geom_center(2) " +" "+ str(letrasx) +" "+ str(letrasy) +" "+ str(letrasz) +" "+ "\n")
                        conta_r = 0

                    if conta_r == 3:

                        num_r3x = float(r3x)
                        num_r3y = float(r3y)
                        num_r3z = float(r3z)
                        gcx3t = (num_r1x + num_r2x + num_r3x)
                        gcx3 = (gcx3t / 3)
                        gcy3t = (num_r1y + num_r2y + num_r3y)
                        gcy3 = (gcy3t / 3)
                        gcz3t = (num_r1z + num_r2z + num_r3z)
                        gcz3 = (gcz3t / 3)

                        if -99.9999 <= gcx3 <= 999.9999:
                            letrasx = gcx3
                        if -999.9999 <= gcx3 <= -100.0000:
                            letrasx = gcx3
                        if -99.9999 <= gcy3 <= 999.9999:
                            letrasy = gcy3
                        if -999.9999 <= gcy3 <= -100.0000:
                            letrasy = gcy3
                        if -99.9999 <= gcz3 <= 999.9999:
                            letrasz = gcz3
                        if -999.9999 <= gcz3 <= -100.0000:
                            letrasz = gcz


                        l_root_f.write("R geom_center(3) " +" "+ str(letrasx) +" "+ str(letrasy) +" "+ str(letrasz) +" "+ "\n")
                        # p_root_f.write(p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                        conta_r = 0
                l_root_f.write(tipo +" "+ l_atomo +" "+ l_nome_atomo + " " + lx + " " + ly + " " + lz + " " + l_enter + "\n")
        if tipo != 'R' and tipo != 'B' and tipo != '\0':
            print("ERRO! |" + tipo + "|")

    l_root.close() 
    l_root_f.close()


# In[7]:


def angulos():
    import math
    bond_count = 0
    #--------variaveis para carregar a proteina
    # raiz
    rp_tipo = ''
    rp_atomo = ''
    rp_nome_atomo =  ''
    rp_aminoacido =  ''
    rp_num_aminoacido =  ''
    c1_rpx =  0.0
    c1_rpy =  0.0
    c1_rpz =  0.0
    rpx =  0.0
    rpy =  0.0
    rpz =  0.0
    num_rpx = 0.0
    num_rpy = 0.0
    num_rpz = 0.0
    rp_enter =  ''

    # aceitador/doador
    bp_tipo =  ''
    bp_atomo =  ''
    bp_nome_atomo =  ''
    bp_aminoacido =  ''
    bp_num_aminoacido =  ''
    c1_bpx =  ''
    c1_bpy =  ''
    c1_bpz =  ''
    bpx =  ''
    bpy =  ''
    bpz =  ''
    num_bpx = 0.0
    num_bpy = 0.0
    num_bpz = 0.0
    bp_enter = ''

    #---------------------------
    #---------variaveis para carregar o ligante
    rl_tipo =  ''
    rl_atomo =  ''
    rl_nome_atomo =  ''
    rl_resto =  ''
    c1_rlx =  ''
    c1_rly =  ''
    c1_rlz =  ''
    rlx =  ''
    rly =  ''
    rlz =  ''
    num_rlx = 0.0
    num_rly = 0.0
    num_rlz = 0.0
    rl_enter =  ''
    bl_tipo =  ''
    bl_atomo =  ''
    bl_nome_atomo =  ''
    bl_resto =  ''
    c1_blx =  ''
    c1_bly =  ''
    c1_blz =  ''
    blx =  ''
    bly =  ''
    blz =  ''
    num_blx = 0.0
    num_bly = 0.0
    num_blz = 0.0
    bl_enter =  ''

    # ----------- saida-arquivos-temporarios com ligante e proteina (atomos mais proximos) -----------
    p_result = open("p_result.tmp", "w")
    l_result = open("l_result.tmp", "w")
    distances = open("dist_result.tmp", "w")

    # -----------------------
    proteina = open("P_ROOT.tmp", "r")
    pligante = open("L_ROOT.tmp", "r")
    
    for p1, p2 in zip(proteina, pligante):       
        p1 = p1.split()
        p2 = p2.split()

        if(p1[1].startswith('geom_center')):
            rp_tipo = p1[0]
            rp_atomo = '-1'
            rpx = p1[2]
            rpy = p1[3]
            rpz = p1[4]
            rp_nome_atomo = ' '
            rp_aminoacido = ' '
            rp_num_aminoacido = ' '
            rp_enter = ' '        
        else:
            rp_tipo = p1[0]
            rp_atomo = p1[1]
            rp_nome_atomo = p1[2]
            rp_aminoacido = p1[3]
            rp_num_aminoacido = p1[4]
            rpx = p1[5]
            rpy = p1[6]
            rpz = p1[7]
            rp_enter = p1[8]

        bp_tipo = rp_tipo
        bp_atomo = rp_atomo
        bp_nome_atomo = rp_nome_atomo
        bp_aminoacido = rp_aminoacido
        bp_num_aminoacido = rp_num_aminoacido
        bpx = rpx
        bpy = rpy
        bpz = rpz
        bp_enter = rp_enter
        
        # 1 N 18.2300 4.9030 9.7160 N.3 1 GLY1 -0.1210
        if(p2[1].startswith('geom_center')):
            rl_tipo = p2[0]
            rl_atomo = '-1'
            rlx = p2[2]
            rly = p2[3]
            rlz = p2[4]
            rl_nome_atomo = ''
            rl_enter = ''
            rl_resto = ''
        else:
            rl_tipo = p2[0]
            rl_atomo = p2[1]
            rl_nome_atomo = p2[2]
            rl_resto = p2[6]
            rlx = p2[3]
            rly = p2[4]
            rlz = p2[5]
            rl_enter = '' 

        bl_tipo = rl_tipo
        bl_atomo = rl_atomo
        bl_nome_atomo = rl_nome_atomo
        bl_resto = rl_resto
        blx = rlx
        bly = rly
        blz = rlz
        bl_enter = rl_enter

        # Copiar strings para outra variável antes de começar a alterar o XYZ
        # Raiz da proteína
        c1_rpx = rpx
        c1_rpy = rpy
        c1_rpz = rpz
        
        # Raiz do ligante
        c1_rlx = rlx
        c1_rly = rly
        c1_rlz = rlz
        
        # Aceitador / doador da proteína
        c1_bpx = bpx
        c1_bpy = bpy
        c1_bpz = bpz
        
        # Aceitador / doador do ligante
        c1_blx = blx
        c1_bly = bly
        c1_blz = blz

        # transforma os 4 vetores em numeros
        # proteina
        # RP = C1
        num_rpx = float(rpx)
        num_rpy = float(rpy)
        num_rpz = float(rpz)
        # ADP = A
        num_bpx = float(bpx)
        num_bpy = float(bpy)
        num_bpz = float(bpz)
        # ligante
        # RL = C2
        num_rlx = float(rlx)
        num_rly = float(rly)
        num_rlz = float(rlz)
        # ADL = B
        num_blx = float(blx)
        num_bly = float(bly)
        num_blz = float(blz)

        # A-B-X
        ABXt1 = num_blx - num_bpx
        ABXt2 = ABXt1 * ABXt1

        # A-B-Y
        ABYt1 = num_bly - num_bpy
        ABYt2 = ABYt1 * ABYt1

        # A-B-Z
        ABZt1 = num_blz - num_bpz
        ABZt2 = ABZt1 * ABZt1

        # soma AB
        ABt1 = ABXt2 + ABYt2
        ABt2 = ABt1 + ABZt2

        # D entre a e b
        DAB = math.sqrt(ABt2)

        # ANGULO-UM
        # A-C1-X
        AC1Xt1 = num_rpx - num_bpx
        AC1Xt2 = AC1Xt1 * AC1Xt1

        # A-C1-Y
        AC1Yt1 = num_rpy - num_bpy
        AC1Yt2 = AC1Yt1 * AC1Yt1

        # A-C1-Z
        AC1Zt1 = num_rpz - num_bpz
        AC1Zt2 = AC1Zt1 * AC1Zt1

        # soma AC1
        AC1t1 = AC1Xt2 + AC1Yt2
        AC1t2 = AC1t1 + AC1Zt2

        # D entre a e C1
        DAC1 = math.sqrt(AC1t2)

        # B-C1-X
        BC1Xt1 = num_rpx - num_blx
        BC1Xt2 = BC1Xt1 * BC1Xt1

        # B-C1-Y
        BC1Yt1 = num_rpy - num_bly
        BC1Yt2 = BC1Yt1 * BC1Yt1

        # B-C1-Z
        BC1Zt1 = num_rpz - num_blz
        BC1Zt2 = BC1Zt1 * BC1Zt1

        # soma BC1
        BC1t1 = BC1Xt2 + BC1Yt2
        BC1t2 = BC1t1 + BC1Zt2

        # D entre a e C1
        DBC1 = math.sqrt(BC1t2)

        # ANGULO-DOIS
        # A-C2-X
        AC2Xt1 = num_rlx - num_bpx
        AC2Xt2 = AC2Xt1 * AC2Xt1

        # A-C2-Y
        AC2Yt1 = num_rly - num_bpy
        AC2Yt2 = AC2Yt1 * AC2Yt1

        # A-C2-Z
        AC2Zt1 = num_rlz - num_bpz
        AC2Zt2 = AC2Zt1 * AC2Zt1

        # soma AC2
        AC2t1 = AC2Xt2 + AC2Yt2
        AC2t2 = AC2t1 + AC2Zt2

        # D entre a e C2
        DAC2 = math.sqrt(AC2t2)

        # B-C2-X
        BC2Xt1 = num_rlx - num_blx
        BC2Xt2 = BC2Xt1 * BC2Xt1

        # B-C2-Y
        BC2Yt1 = num_rly - num_bly
        BC2Yt2 = BC2Yt1 * BC2Yt1

        # B-C2-Z
        BC2Zt1 = num_rlz - num_blz
        BC2Zt2 = BC2Zt1 * BC2Zt1

        # soma BC2
        BC2t1 = BC2Xt2 + BC2Yt2
        BC2t2 = BC2t1 + BC2Zt2

        # D entre a e C2
        DBC2 = math.sqrt(BC2t2)

        # ANGULOS
        DAB_2 = DAB * DAB

        # ANGULO-UM
        DAC1_2 = DAC1 * DAC1
        DBC1_2 = DBC1 * DBC1

        # soma AB AC1
        S_DAB_DAC1_DBC1_t1 = DAB_2 + DAC1_2
        S_DAB_DAC1_DBC1_t2 = S_DAB_DAC1_DBC1_t1 - DBC1_2

        # AB * AC1
        r2V_DAB_DAC1_t1 = DAB * DAC1
        r2V_DAB_DAC1_t2 = 2 * r2V_DAB_DAC1_t1
        if r2V_DAB_DAC1_t2 == 0:
            r2V_DAB_DAC1_t2 = 0.00001
        # angulo 1
        r2V_DAB_DAC1_t3 = S_DAB_DAC1_DBC1_t2 / r2V_DAB_DAC1_t2
        angulo_1_t1 = math.acos(r2V_DAB_DAC1_t3)

        # radianos em graus
        angulo_1_t2 = angulo_1_t1 * 180
        angulo_1 = angulo_1_t2 / 3.141618

        # ANGULO-DOIS
        DAC2_2 = DAC2 * DAC2
        DBC2_2 = DBC2 * DBC2

        # soma AB AC2
        S_DAB_DAC2_DBC2_t1 = DAB_2 + DBC2_2
        S_DAB_DAC2_DBC2_t2 = S_DAB_DAC2_DBC2_t1 - DAC2_2

        # AB * AC2
        r2V_DAB_DAC2_t1 = DAB * DBC2
        r2V_DAB_DAC2_t2 = 2 * r2V_DAB_DAC2_t1

        # angulo 2
        if r2V_DAB_DAC2_t2 == 0:
            r2V_DAB_DAC2_t2 = 0.0000001
        
        r2V_DAB_DAC2_t3 = S_DAB_DAC2_DBC2_t2 / r2V_DAB_DAC2_t2
        angulo_2_t1 = math.acos(r2V_DAB_DAC2_t3)

        # radianos em graus
        angulo_2_t2 = angulo_2_t1 * 180
        angulo_2 = angulo_2_t2 / 3.141618

        # SALVA ARQUIVO PDB COM ATOMOS E RAIZES OS ATOMOS MAIS PROXIMOS
        num1_bl_atomo = 0
        num1_bp_atomo = 0
        num_bl_atomo = 0
        num_bp_atomo = 0

        num_bl_atomo = int(bl_atomo)
        num_bp_atomo = int(bp_atomo)

        if angulo_1 >= 60 and angulo_2 >= 60:
            # contador de ligações de hidrogênio
            bond_count += 1
            
            # -------------- saída ligante TEMPORÁRIA
            num1_bl_atomo = num_bl_atomo
            l_result.write("HETATM" + " " + bl_atomo + " " + bl_nome_atomo + " BLK " + bl_resto + " " + c1_blx + " " + c1_bly + " " + c1_blz + " " + bl_enter+'\n')

            # ------------------------------------------
            # -------------- saída proteína TEMPORÁRIA
            num1_bp_atomo = num_bp_atomo
            p_result.write("ATOM" + " " + bp_atomo + " " + bp_nome_atomo + " " + bp_aminoacido + " " + bp_num_aminoacido + " " + c1_bpx + " " + c1_bpy + " " + c1_bpz + " " + bp_enter+'\n')
            # ------------------------------------------
            distance = DAB
            distances.write(str(distance) + " " + "\n")

    p_result.close()  # fecha a saída da proteína
    l_result.close()  # fecha a saída do ligante

    distances.close()

    # mostra a contagem total de pontes de hidrogênio
    """
    bound_total = open("bond_count.txt", "w")
    bound_total.write(str(bond_count))
    bound_total.close()
    """

    # -------------- saída: resultados detalhados
    # entradas da proteína e ligante
    p_line = ""
    l_line = ""
    b_distances = ""

    # gera arquivo de log
    bonds_log = open("bonds.log", "w")



    bonds_log.write("Hydrogen bonds found: " + str(bond_count) + "\n\n")
    hydrogen_B = bond_count

    p_bonds = open("p_result.tmp", "r")
    l_bonds = open("l_result.tmp", "r")
    

    with open("p_result.tmp", "r") as p_bonds, open("l_result.tmp", "r") as l_bonds:
        for proteina, ligante in zip(p_bonds, l_bonds):
            p_line = p_bonds.readline().strip()
            l_line = l_bonds.readline().strip()
            # b_distances = d_bonds.readline().strip()
            bonds_log.write(p_line + l_line)
            # bonds_log.write("distance: " + b_distances)
            bonds_log.write("\n")

    p_bonds.close()
    l_bonds.close()
    
    bonds_log.write("--------------------------------------------------------\n")
    bonds_log.write("        TABLE: INTERMOLECULAR HYDROGEN BONDS            \n")
    bonds_log.write("\n")
    bonds_log.write("            Protein Ligand Distance(A)\n                  ")

    d_bonds = open("dist_result.tmp", "r")
    p_bonds_final = open("p_result.tmp", "r")
    l_bonds_final = open("l_result.tmp", "r")

    for p_line in p_bonds_final:
        l_line = l_bonds_final.readline()
        distance = d_bonds.readline().strip()
        bonds_log.write(p_line.strip() + "  " + l_line.strip() + "  " + distance + "\n")
    
    # Reset the read position to the beginning of the file
    p_bonds_final.close()
    p_bonds_final.close()
    d_bonds.close()
    
    d_bonds = open("dist_result.tmp", "r")
    p_bonds_final = open("p_result.tmp", "r")
    l_bonds_final = open("l_result.tmp", "r")
    # abre arquivo para depois comparar o limite máximo de ligações de hidrogênio
    limit_l = open("limit_l.tmp", "w")
    limit_p = open("limit_p.tmp", "w")
    
    for p_line, l_line, b_distances in zip(p_bonds_final, l_bonds_final, d_bonds):
        p_line = p_line.strip().split(" ")
        l_line = l_line.strip().split(" ")
        if str(p_line[3]):
            aux = ( str(p_line[3]) +" "+str(p_line[4]) +" "+str(p_line[8]) + " " + str(p_line[2]) + "\n")
            bonds_log.write(aux)
            limit_p.write( aux )

        # Proteína: divisor
        bonds_log.write(" ")

        aux = ( str(l_line[3]) +" "+str(l_line[4]) +" "+str(l_line[1]) + "\n")
        bonds_log.write(aux)
        limit_l.write( aux )

        # Ligante: divisor
        bonds_log.write("         ")

        # Distância
        bonds_log.write(b_distances)
        bonds_log.write("\n")

    p_bonds_final.close()
    l_bonds_final.close()
    d_bonds.close()
    limit_p.close()
    limit_l.close()
    bonds_log.close()
    
    return bond_count


# In[8]:


def raiz_ligante(protein_bond,protein_tripos):
    # print (protein_tripos)
    PL_DIST = {}
    f = open('PL_DIST.tmp','r')   
    for l in f:
        l = l.strip().split(' ')
        PL_DIST[l[0]] = l
    f.close()
    outlt = open('limit_type.tmp','w')
    out = open('L_ROOT_0.tmp','w')
    for bl in (protein_bond):
        if len(bl) > 2:
            if bl[1] in PL_DIST:
                out.write('R ' + str(' '.join(protein_tripos[bl[2]]) + '\n'))
                out.write('B ' + str(' '.join(PL_DIST[bl[1]]) + '\n'))
                outlt.write(str(protein_tripos[bl[2]][0]) + '\n')
            if bl[2] in PL_DIST:
                out.write('R ' + str(' '.join(protein_tripos[bl[1]]) + '\n'))
                out.write('B ' + str(' '.join(PL_DIST[bl[2]]) + '\n'))
                outlt.write(str(protein_tripos[bl[2]][0]) + '\n')
    out.close()
    outlt.close()


# In[9]:


def saida_PDB(protein):

    # Saída final
    v = 'file'
    hb_file_name = v + "_H-Bonds.pdb"
    result_PDB = open(hb_file_name, "w")

    # Arquivos de entrada (TEMP)
    p_result = open("p_result.tmp", "r")
    l_result = open("l_result.tmp", "r")

    # Variáveis de leitura
    ATOM = ""
    esp1 = ""
    esp2 = ""
    esp3 = ""
    bp_enter = ""
    bp_atomo = ""
    bp_nome_atomo = ""
    bp_aminoacido = ""
    bp_num_aminoacido = 0
    num_bp_num_aminoacido = 0
    num_bp_atomo = 0
    atomo_anterior = 0
    num_p_num_aminoacido = 0
    num_p_atomo = 0
    bpx = ""
    bpy = ""
    bpz = ""
    p_atomo = ""
    p_nome_atomo = ""
    p_aminoacido = ""
    p_num_aminoacido = ""
    px = ""
    py = ""
    pz = ""
    p_enter = ""
    num_aminoacido_anterior = 0

    while True:
        line = p_result.readline()
        if not line:
            break
        # ['ATOM', '-1', '', '', '', '17.092', '-0.025499999999999995', '8.882000000000001']
        l = line.strip().split(' ')

        if l[1] == '-1':
            # Faz a leitura dos campos
            ATOM = l[0]
            bp_atomo = l[1]
            esp1 = ''
            bp_nome_atomo = ''
            bp_aminoacido = ''
            bp_num_aminoacido = 0
            esp2 = ''
            bpx = l[5]
            bpy = l[6]
            bpz = l[7]
            esp3 = ''
            bp_enter = ''
            num_bp_num_aminoacido = int(bp_num_aminoacido)
            num_bp_atomo = int(bp_atomo)
        else:
            # ATOM 81 OG SER 12 10.019 -0.09 20.273 A
            # Faz a leitura dos campos
            ATOM = l[0]
            bp_atomo = l[1]
            esp1 = ''
            bp_nome_atomo = l[2]
            bp_aminoacido = l[3]
            bp_num_aminoacido = l[4]
            esp2 = ''
            bpx = l[5]
            bpy = l[6]
            bpz = l[7]
            esp3 = ''
            bp_enter = ''

            num_bp_num_aminoacido = int(bp_num_aminoacido)
            num_bp_atomo = int(bp_atomo)

        for p_line in protein:
            # Faz a leitura dos campos
            #   1    2     3     4    5       6         7       8        9       10     11
            # ['1', 'N', 'GLY', 'A', '1', '18.230', '4.903', '9.716', '1.00', '28.23', 'N']
            # ['2', 'CA', 'GLY', 'A', '1', '16.795', '5.303', '9.674', '1.00', '26.95', 'C']
            p_atomo = p_line[0]
            p_nome_atomo = p_line[1]
            p_aminoacido = p_line[2]
            p_num_aminoacido = p_line[4]
            px = p_line[5]
            py = p_line[6]
            pz = p_line[7]
            p_enter = p_line[8]

            num_p_num_aminoacido = int(p_num_aminoacido)
            num_p_atomo = int(p_atomo)

            if num_p_num_aminoacido == num_bp_num_aminoacido:
                if num_p_atomo < num_bp_atomo:
                    result_PDB.write("ATOM" +" "+ p_atomo +" "+ esp1 +" "+ p_nome_atomo +" "+ p_aminoacido +" "+ p_num_aminoacido +" "+ esp2 +" "+ px +" "+ py +" "+ pz +" "+ esp3 +" "+ p_enter + "\n")
                elif num_p_atomo == num_bp_atomo:
                    result_PDB.write("ATOM" +" "+ p_atomo +" "+ esp1 +" "+ p_nome_atomo +" "+ p_aminoacido +" "+ p_num_aminoacido +" "+ esp2 +" "+ px +" "+ py +" "+ pz +" "+ esp3 +" "+ p_enter + "\n")
                elif num_p_atomo > num_bp_atomo:
                    result_PDB.write("ATOM" +" "+ p_atomo +" "+ esp1 +" "+ p_nome_atomo +" "+ p_aminoacido +" "+ p_num_aminoacido +" "+ esp2 +" "+ px +" "+ py +" "+ pz +" "+ esp3 +" "+ p_enter + "\n")

        num_aminoacido_anterior = num_bp_num_aminoacido

    while True:
        ch = l_result.read(1)
        if not ch:
            break
        result_PDB.write(ch)

    # Fecha os arquivos
    result_PDB.close()
    l_result.close()
    p_result.close()


# In[10]:


def salva_proteina():
    L = 0
    lines_total = 0
    v = []
    p_result = open("p_result.tmp", "r")  # Substitua "p_result.txt" pelo nome do arquivo apropriado
    v = []
    v.clear()
    ligand_name = "l_result.tmp"  # Substitua pelo nome adequado
    #v = ligand_name.split('')
    hb_file_name = "file_H-Bonds.pdb"
    
    p_result = open(hb_file_name, "r")
    
    lines = p_result.readlines()
    lines_total = len(lines)
    p_result.close()
    
    vec = [""] * lines_total
    
    p_result2 = open(hb_file_name, "r")
    
    for line in p_result2:
        vec[L] = line
        L += 1
    p_result2.close()
    
    #del_rep(vec)
    
    p_result3 = open(hb_file_name, "w")
    for i in range(lines_total):
        p_result3.write(vec[i])
    p_result3.close()


# In[23]:


def calcula_RT(ligand):
    limit_type = open("limit_type.tmp",'r')
    line_count = 0
    current_line = 0
    atom_1 ={}
    
    for line in limit_type:
        line_count += 1
        atom_1[current_line] = line.strip()
        current_line += 1

    with open("RT_found.tmp") as RT_found:
        rt_string = []
        for line in RT_found:
            rt_string.append(line.strip())

    limit_type.close()
    
    # Verificar
    lig_line = 0
    RT_found = open("RT_found.tmp")
    rt_line = {}
    for line in RT_found:
        if not line:
            break
        rt_line[len(line)] = line
        lig_line += 1
    total_lig_line = lig_line
    rt_string = {}
    lig_line = 0
    RT_found.close()
    
    RT_found = open("RT_found.tmp")
    for line in RT_found:
        if not line:
            break
        rt_string[lig_line] = line.strip()
        lig_line += 1
    i = 0
    j = 0
    aux = rt_string.copy()
    for i in range(lig_line):
        for j in range(lig_line):
            # Entender essa matrix
            if i != j and rt_string[i] == aux[j]:
                aux[i] = 'REP'
    lines_total2 = 0
    lines_total2 -= 1
    j = 0
    # I'm here ()
    #rt_string_copy = bytearray(total_lig_line*8)
    rt_string_copy = list()
    for i in range(lig_line):
        if aux[i] != 'REP':
            indice = int(j)
            rt_string_copy.append(aux[i])
            j += 1
            lines_total2 += 1
    current_line_0 = 0
    current_line_1 = 0
    lig_line = 0
    
    type_l = 0
    mol2_lines = list()
    for line in ligand:
        mol2_lines.append(line)
        type_l += 1

    mol2_type = {}
    
    mol2_type = list()
    type_l = 0
    for mol2_line in mol2_lines:
        mol2_type.append(mol2_line)
        type_l += 1

    sp2_count = 0
    atom_count = 0
    marker = 0
    pos = 0
    pair_count1 = 0
    pair_count2 = 1
    num_atom1 = 0
    num_atom2 = 0
    atom1 = ''
    atom2 = ''
    number_count = 0
    pair = 0
    next = 0
    s_bond = 0
    RT = 0.0
    # Verificar limit_type.tmp file
    while current_line_0 < lines_total2:
        if marker < line_count and current_line_0 < lines_total2:
            while marker < line_count:
                number_count += 1
                next = marker + 1

                if rt_string[current_line_0] == atom_1[marker]:
                    atom_count += 1
                    try:
                        if atom_1[next][0] == '1':
                            atom1 = atom_1[marker]
                            num_atom1 = int(atom1)
                            if mol2_type[num_atom1][0] == '2':
                                sp2_count += 1
                            if number_count % 2 != 0:
                                pair = 2
                            if number_count % 2 == 0:
                                pair = -2
                            atom2 = atom_1[marker + pair]
                            num_atom2 = int(atom2)
                            if mol2_type[num_atom2][0] == '2':
                                sp2_count += 1
                            if sp2_count < 2 and mol2_type[num_atom2][7] != 'H' and mol2_type[num_atom1][7] != 'H':
                                s_bond += 1
                    except:
                        marker += 2
                        break
                marker += 2
        marker = 0
        
        if atom_count > 1:
            if s_bond == 0:
                RT += 0
            elif s_bond == 1 or s_bond >= 3:
                RT += 0.5
            elif s_bond == 2:
                RT += 1

        s_bond = 0
        marker = 0
        sp2_count = 0
        atom_count = 0
        current_line_0 +=1
        
    return RT


# In[12]:


def limit():
    p_list = open("limit_p.tmp", "r")
    line_p_count = 0
    line_l_count = 0
    line_d_count = 0
    p_position = 0
    l_position = 0
    d_position = 0
    d_smaller = 0
    copied = 0
    compared = 0
    MAX_O = 4
    MAX_C = 4
    MAX_P = 4
    MAX_S = 4
    MAX_F = 4
    MAX_N = 4
    MAX = 0
    ch = ''
    TOTAL = 0
    p_TOTAL = 0

    for line in p_list:
        ch = line.split(' ')[0]
        #if ch[21] == '\n':
        p_TOTAL += 1

    p_list.close()
    p_list = open("limit_p.tmp", "r")
    
    ch_p = {}
    
    line_p = 0
    for line in p_list:
        ch = line.split(' ')[0]
        ch_p[line_p] = ch
        line_p += 1
    line_p = 0
    line_p_count = 0
    p_list.close()
    p_root = open("P_ROOT_0.tmp", "r")
    p_liged = 0
    #p_liged_list = [[''] * 3 for _ in range(p_TOTAL)]
    p_liged_list = {}
    p_file_line = 0
    for pfl in p_root:
        # ['B', '19', 'N', 'ASP', '4', '17.098', '0.654', '8.646', 'A\n']
        p_root_line = pfl.strip().split(" ")
        if p_root_line[0] == 'R' and p_root_line[8] != 'H':
            p_liged += 1
        if p_root_line[0] == 'B':
            p_file_line += 1
        if (p_root_line[0] == 'B' and ( p_file_line - 1 != line_p_count or p_root_line[0] != ch_p[line_p_count]) ):
            p_liged = 0
        
        p_liged_list[line_p_count] = str(p_liged)
        p_liged = 0
        line_p_count += 1
        p_file_line = 0
    p_root.close()
    
    #while line_p_count < p_TOTAL:
    #    p_function()
    
    p_root.close()
    l_TOTAL = 0
    l_list = open("limit_l.tmp", "r")
    for line in l_list:
        l_TOTAL += 1
        # BLK N.am 12
        
    l_list.close()
    
    l_list = open("limit_l.tmp", "r")
    ch_l = [[''] * 17 for _ in range(l_TOTAL)]    
    line_l = 0
    
    for line in l_list:
        ch_l[line_l] = line.strip().split(' ')[1]
        line_l += 1

    line_l_count = 0  # a linha que esta sendo comparada com as outras

    l_list.close()
    
    l_root = open("L_ROOT_0.tmp")
    l_file_line = 0
    l_liged = 0
    l_liged_list = {}
    l_root_line = list()
    # --------------------------------------------
    # l_function:
    # --------------------------------------------
    #while not l_root.eof() and (l_file_line - 1) < line_l_count:
    for lr in l_root:
        l_root_line = lr.strip().split()
        if l_root_line[0] == 'R' and l_root_line[6] != 'H' and l_root_line[6] != 'h':
        #if l_root_line[0] == 'R':
            l_liged += 1
            l_liged_list[line_l_count] = l_liged
        if l_root_line[0] == 'B':
            l_file_line += 1
            l_liged_list[line_l_count] = l_liged
        if (l_root_line[0] == 'B' and ((l_file_line - 1) != line_l_count or l_root_line[0] != ch_l[line_l_count])):
            l_liged = 0
            l_liged_list[line_l_count] = l_liged
        line_l_count += 1

    l_liged = 0
    
    l_file_line = 0
    l_root.close()
    # -------------------------abre distancias-------------------
    dist_line = {}
    ch3 = [None] * 12
    line_d = 0
    d_list = open("dist_result.tmp")
    #while not d_list.eof() and line_d < p_TOTAL:
    for dl in d_list:
        ch3 = dl.strip().split(' ')
        dist_line[line_d] = ch3
        line_d += 1

    line_d = 0
    line_d_count = 0

    liged = 0  # numero de atomos ja ligados
    found = 0
    string_line_count = 0
    limit_type = open("limit_type.tmp")
    type = [None] * 8
    atom = [None] * 8
    refined = 0
    while line_p_count < p_TOTAL:
        liged = int(p_liged_list[line_p_count])
        """
        if ch_p[line_p_count][0] == 'C':
            MAX = MAX_C - liged  # cout << endl << "C" << endl;
        if ch_p[line_p_count][0] == 'N':
            MAX = MAX_N - liged  # cout << endl << "N" << endl;
        if ch_p[line_p_count][0] == 'P':
            MAX = MAX_P - liged  # cout << endl << "P" << endl;
        if ch_p[line_p_count][0] == 'O':
            MAX = MAX_O - liged  # cout << endl << "O" << endl;
        if ch_p[line_p_count][0] == 'S':
            MAX = MAX_S - liged  # cout << endl << "S" << endl;"""

        liged = int(p_liged_list[line_p_count])
        if ch_p[line_p_count][0].startswith('C'):
            MAX = MAX_C - liged  # cout << endl << "C" << endl;
        if ch_p[line_p_count][0].startswith('N'):
            MAX = MAX_N - liged  # cout << endl << "N" << endl;
        if ch_p[line_p_count][0].startswith('P'):
            MAX = MAX_P - liged  # cout << endl << "P" << endl;
        if ch_p[line_p_count][0].startswith('O'):
            MAX = MAX_O - liged  # cout << endl << "O" << endl;
        if ch_p[line_p_count][0].startswith('S'):
            MAX = MAX_S - liged  # cout << endl << "S" << endl;
            
        while string_line_count < p_TOTAL:
            if (
                ch_p[line_p_count][12] != 'X'
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
                and ch_p[line_p_count] == ch_p[string_line_count]
            ):
                found += 1  # encontrou um atomo igual ao da string selecionada
                if found > MAX:
                    ch_p[string_line_count] = 'X'  # marca este atomo para excluir da procura
                    ch_l[string_line_count] = 'X'
                    dist_line[string_line_count][0] = 'X'
            string_line_count += 1
        # cout << endl << "F" << found << endl;
        string_line_count = 0  # reseta as linhas que vao serem comparadas com a selecionada
        found = 0  # reseta numero de encontrados
        line_p_count += 1  # seleciona a proxima string a ser comparada com as outras
        # d_list.close()


    line_p_count = 0
    line_l_count = 0
    # -------------------------------------------------------------
    # ------refinamento ligante-------------------------------------------------------
    # ------------------------------------------------------------------
    MAX = 0
    
    for line_l_count in l_liged_list:
        liged = l_liged_list[line_l_count]

        try:
            if ch_l[line_l_count].startswith('C'):
                MAX = MAX_C - liged  # cout << endl << "C" << endl;
            if ch_l[line_l_count].startswith('N'):
                MAX = MAX_N - liged  # cout << endl << "N" << endl;
            if ch_l[line_l_count].startswith('P'):
                MAX = MAX_P - liged  # cout << endl << "P" << endl;
            if ch_l[line_l_count].startswith('O'):
                MAX = MAX_O - liged  # cout << endl << "O" << endl;
            if ch_l[line_l_count].startswith('S'):
                MAX = MAX_S - liged  # cout << endl << "S" << endl;
        except:
            continue
            #print ( f"l_liged_list: {len(l_liged_list)} ch_l: {len(ch_l)}" )
            
        for l in limit_type:
            atom = l.strip()
            type = l.strip()
            # print (f"----------")
            # print (f"Type: {type}")
            # print (f"atom: {atom}")
            # print (f"ch_l: {ch_l[line_l_count]}")
            # print (f"----------")
            # cout << type << "*" << endl;
            try:
                if ( type == '2' and atom[1] == ch_l[line_l_count] ):
                    MAX -= 1
            except:
                MAX -= 0

        limit_type.seek(0)  # move to the start of the file
        string_line_count = 0

        while string_line_count < l_TOTAL:
            if line_l_count < l_TOTAL:
                if (ch_l[line_l_count] != 'X' and ch_l[line_l_count] == ch_l[string_line_count]):
                        found += 1  # encontrou um atomo igual ao da string selecionada
                if found > MAX:
                            ch_p[string_line_count] = 'X'  # marca este atomo para excluir da procura
                            ch_l[string_line_count] = 'X'
                            dist_line[string_line_count] = 'X'
            string_line_count += 1

        # cout << endl << "F" << found << "*" << ch_l[line_l_count][9] << endl;
        string_line_count = 0  # reseta as linhas que vao serem comparadas com a selecionada
        found = 0  # reseta numero de encontrados
        line_l_count += 1  # seleciona a proxima string a ser comparada com as outras
        # d_list.close()

    line_p_count = 0
    line_l_count = 0
    string_line_count = 0
    while string_line_count < l_TOTAL:
        if ch_l[string_line_count] != 'X':
            refined += 1
        string_line_count += 1

    string_line_count = 0


# In[13]:


def result_score_calc(protein_name,ligand_pdb,aux,seqaa,hydrogen_B,aa_table):

    import os
    import subprocess
    HC_total2 = aux[0]
    VDW_total = aux[1]
    RT = aux[2]
    #hydrogen_B = aux[3]
    repulsive = aux[3]
    london = aux[4]
    # HC_total2, VDW_total, RT, hydrogen_B, repulsive, london
    global ligand_name, ligand_type
    #global hydrophobicity, surface_tension

    HB = 0
    HB = refined
    
    result_file_name = ligand_name.split('.')[0] + "_result.txt"

    with open(result_file_name, 'a') as result_score:
        result_score.write(f"Hydrophobic contacts:        {HC_total2}\n")
        result_score.write(f"Van der waals:               {VDW_total}\n")
        result_score.write(f"Deformation effect:          {RT}\n")
        result_score.write(f"Hydrogen bonds (HB):         {hydrogen_B}\n")
        result_score.write(f"Repulsive VDW score:         {repulsive}\n")
        result_score.write(f"London dispersion force:     {london}\n")
    
    test  = 1
    test2 = 5
    result = summation(test, test2)

    #ASA1, ASA2 = 0, 0
    ASA1 = (subprocess.run(["python", "asa/asa.py", protein_name], stdout=subprocess.PIPE).stdout)
    ASA2 = (subprocess.run(["python", "asa/asa.py", ligand_pdb], stdout=subprocess.PIPE).stdout)

    import re
    ASA1 = re.search(r'([\d.]+)', ASA1.decode('utf-8'))
    ASA1 = float(ASA1.group(0))
      
    ASA2 = re.search(r'([\d.]+)', ASA2.decode('utf-8'))
    ASA2 = float(ASA2.group(0))

    hydro_file = "hydrophobicity.param"
    tension_file = "tension.param"

    # Cálculo da hidrofobicidade e da tensão de superfície
    hydro_map = create_map_tables(hydro_file)
    tension_map = create_map_tables(tension_file)

    protein_file_name = "protein.tmp"
    p_dist_file_name = "PP_DIST.tmp"

    total_hydrophobicity = 0
    total_surface_tension = 0
    contact_hydrophobicity = 0
    contact_surface_tension = 0
    for a in seqaa:
        total_hydrophobicity += hydro_map[aa_table[a]]
        total_surface_tension += tension_map[aa_table[a]]

    contact_hydrophobicity,contact_surface_tension = surface_tension_hydrophobicity_calculator(hydro_map, tension_map, p_dist_file_name)
   
    # New features
    r1 = calculate_proportion(seqaa)
    r2 = CalculateAAComposition(seqaa)
    r3 = CalculateDipeptide(seqaa)
    r4 = CalculateTripeptide(seqaa)

    res = list()
    
    #out = open(str3, 'w')
    key_ = ['HC_total2','VDW_total','RT','hydrogen_B','ASA1','ASA2','repulsive','london','contact_hydrophobicity','total_hydrophobicity','contact_surface_tension','total_surface_tension']

    res.append(HC_total2)
    res.append(VDW_total)
    res.append(RT)
    res.append(hydrogen_B)
    res.append(ASA1)
    res.append(ASA2)
    res.append(repulsive)
    res.append(london)
    res.append(contact_hydrophobicity)
    res.append(total_hydrophobicity)
    res.append(contact_surface_tension)
    res.append(total_surface_tension)
    for f in r1:
        key_.append(f)
        res.append(r1[f])
    for f in r2:
        key_.append(f)
        res.append(r2[f])
    for f in r3:
        key_.append(f)
        res.append(r3[f])
    for f in r4:
        key_.append(f)
        res.append(r4[f])

    f = open('features.txt','w')
    f.write('\t'.join(map(str, key_)))
    f.close()
    
    return res,key_


# In[14]:


def create_map_tables(param_file_name):
    param_map = {}
    
    with open(param_file_name, 'r') as param_file:
        for param_line in param_file:
            param_v = param_line.strip().split(':')
            param_map[param_v[0]] = float(param_v[1])
    
    return param_map


# In[15]:


def summation(a, b):
    x = a
    _sum = 0

    while x <= b:
        _sum += x
        x += 1

    return _sum


# In[16]:


def surface_tension_hydrophobicity_calculator(hydro_map, tension_map, infile_name):
    aa_numbers = []
    counter = 0
    hydrophobicity = 0
    surface_tension = 0

    with open(infile_name, 'r') as infile:
        for p_dist_line in infile:

            v = p_dist_line.strip().split(' ')[2]
            if hydro_map[v]:
                hydrophobicity += hydro_map[v]
            if tension_map[v]:
                surface_tension += tension_map[v]

    return hydrophobicity, surface_tension


# In[17]:


def CalculateTripeptide(aaseq):
    results = {}
    for aa in AATripeptide:
        results[aa] = 0 if aaseq.count(aa) == 0 else round(float(aaseq.count(aa)) / float(len(aaseq)), 5)
    return results


# In[18]:


def CalculateDipeptide(aaseq):
    results = {}
    for i in AA:
        for j in AA:
            Dipeptide = i+j
            results[Dipeptide] = 0 if aaseq.count(Dipeptide) == 0 else round(float(aaseq.count(Dipeptide)) / float(len(aaseq)), 5)            
    return results


# In[19]:


def CalculateAAComposition(aaseq):
    results = {}
    for i in AA:
        results[i] = 0 if aaseq.count(i) == 0 else round(float(aaseq.count(i)) / float(len(aaseq)), 5)
    return results


# In[20]:


def calculate_proportion(seq):
    # **************************** Take a important note *******************************
    # Selenocysteine is the Se-analogue of cysteine. (U)
    # Pyrrolysine is synthesized in vivo by joining two molecules of L-lysine. (O)
    # **********************************************************************************
    mapVolume = {}
    mapMass = {}
    mapHydro = {}
    
    mapVolume['A'] = 88.6
    mapVolume['R'] = 173.4
    mapVolume['N'] = 114.1
    mapVolume['D'] = 111.1
    mapVolume['C'] = 108.5
    mapVolume['Q'] = 143.8
    mapVolume['E'] = 138.4
    mapVolume['G'] = 60.1
    mapVolume['H'] = 153.2
    mapVolume['I'] = 166.7
    mapVolume['L'] = 166.7
    mapVolume['K'] = 168.6
    mapVolume['M'] = 162.9
    mapVolume['F'] = 189.9
    mapVolume['P'] = 112.7
    mapVolume['S'] = 89.0
    mapVolume['T'] = 116.1
    mapVolume['W'] = 227.8
    mapVolume['Y'] = 193.6
    mapVolume['V'] = 140.0
    mapVolume['O'] = 2 * 166.7
    mapVolume['U'] = 108.5 # The same value that C
    mapVolume['B'] = 111.1
    mapVolume['Z'] = 138.4
    mapVolume['J'] = 162.9

    mapMass['A'] = 89.0
    mapMass['R'] = 174.0
    mapMass['N'] = 132.0
    mapMass['D'] = 133.0
    mapMass['C'] = 121.0
    mapMass['Q'] = 146.0
    mapMass['E'] = 147.0
    mapMass['G'] = 75.0
    mapMass['H'] = 155.0
    mapMass['I'] = 131.0
    mapMass['L'] = 131.0
    mapMass['K'] = 146.0
    mapMass['M'] = 149.0
    mapMass['F'] = 165.0
    mapMass['P'] = 115.0
    mapMass['S'] = 105.0
    mapMass['T'] = 119.0
    mapMass['W'] = 204.0
    mapMass['Y'] = 181.0
    mapMass['V'] = 117.0
    mapMass['O'] = 237.14773
    mapMass['U'] = 150.95364
    mapMass['B'] = 133.0
    mapMass['Z'] = 147.0
    mapMass['J'] = 196.106

    mapHydro['A'] = 1.8
    mapHydro['R'] = -4.5
    mapHydro['N'] = -3.5
    mapHydro['D'] = -3.5
    mapHydro['C'] = 2.5
    mapHydro['Q'] = -3.5
    mapHydro['E'] = -3.5
    mapHydro['G'] = -0.4
    mapHydro['H'] = -3.2
    mapHydro['I'] = 4.5
    mapHydro['L'] = 3.8
    mapHydro['K'] = -3.9
    mapHydro['M'] = 1.9
    mapHydro['F'] = 2.8
    mapHydro['P'] = -1.6
    mapHydro['S'] = -0.8
    mapHydro['T'] = -0.7
    mapHydro['W'] = -0.9
    mapHydro['Y'] = -1.3
    mapHydro['V'] = 4.2 
    mapHydro['O'] = 2 * 3.8
    mapHydro['U'] = 2.5 # The same value that C
    mapHydro['B'] = -3.5
    mapHydro['Z'] = -3.5
    mapHydro['J'] = 1.9
    
    length = len(seq)
    
    Nonpolar_Aliphatic = 0.0
    Aromatic = 0.0
    Polar_Uncharged = 0.0
    Positively_Charged = 0.0
    Negatively_Charged = 0.0
    mass = 0.0
    volume = 0.0
    hydro = 0.0

    for aa in seq:
        if ( 'G' in aa or 'A' in aa  or 'P'  in aa  or 'V'  in aa  or 'L'  in aa  or 'I'  in aa or 'M'  in aa or 'O'  in aa or 'J'  in aa):
            Nonpolar_Aliphatic += 1
        elif ( 'F'  in aa or 'Y' in aa or 'W'  in aa ):
            Aromatic += 1
        elif ( 'S'  in aa or 'T'  in aa  or 'C'  in aa  or 'N' in aa  or 'Q'  in aa or 'U'  in aa ):
            Polar_Uncharged +=1
        elif ( 'K' in aa or 'H'  in aa or 'R'  in aa ):
            Positively_Charged +=1
        elif ( 'D' in aa or  'E' in aa or 'B' in aa or 'Z' in aa):
            Negatively_Charged +=1
        else:
            print ('')

        if (aa in mapHydro):
            hydro += mapHydro[aa]
            mass += mapMass[aa]
            volume += mapVolume[aa]
    
    Nonpolar_Aliphatic  = 0 if Nonpolar_Aliphatic == 0 else Nonpolar_Aliphatic / len(seq)
    Aromatic = 0 if Aromatic  == 0 else Aromatic / len(seq)
    Polar_Uncharged = 0 if Polar_Uncharged == 0 else Polar_Uncharged / len(seq)
    Positively_Charged  = 0 if Positively_Charged == 0 else Positively_Charged / len(seq)
    Negatively_Charged  = 0 if Negatively_Charged == 0  else Negatively_Charged / len(seq)
    mass= 0 if mass == 0 else mass / len(seq)
    volume  = 0 if volume == 0 else volume / len(seq)
    hydro  = 0 if hydro == 0 else hydro / len(seq)
    
    res = {}
    res['Nonpolar_Aliphatic'] = Nonpolar_Aliphatic
    res['Aromatic'] = Aromatic
    res['Polar_Uncharged'] = Polar_Uncharged
    res['Positively_Charged'] = Positively_Charged
    res['Negatively_Charged'] = Negatively_Charged
    res['mass'] = mass
    res['volume'] = volume
    res['hydro'] = hydro
    
    return res


# In[21]:


def deleta_temp():
    import os
    os.remove("dist_result.tmp")
    os.remove("limit_l.tmp")
    os.remove("limit_p.tmp")
    os.remove("limit_type.tmp")
    os.remove("l_result.tmp")
    os.remove("L_ROOT_0.tmp")
    os.remove("L_ROOT.tmp")
    os.remove("PL_DIST.tmp")
    os.remove("PP_DIST.tmp")
    os.remove("p_result.tmp")
    os.remove("P_ROOT_0.tmp")
    os.remove("P_ROOT.tmp")
    os.remove("RT_found.tmp")    


# In[22]:


import sys

"""
First version.
Author: Jose Cleydson F Silva
Date: 09/29/2323
Convert: jupyter nbconvert --to script  FextractOmics.ipynb
*******************************************************************************************************************************************
New Version;
12/18/2023 - New features and Hidrogen_B feature.
*******************************************************************************************************************************************
01/11/2024: (Raquel and Layla) We are debugging your version of the code and it seems like some of the residue counts don't match:
For example, we expected to see SF having a frequency of 2, which divided by the length of the protein should be 0.00481. Is that right? 
Or how do you calculate those frequencies?
And TA doesn't show up at all.
Fixed: 01/11/24
******************************************************************************************************************************************
# https://www.sciencefriday.com/wp-content/uploads/2018/07/amino-acid-abbreviation-chart.pdf
# https://en.wikipedia.org/wiki/Amino_acid

# Non-proteinogenic and proteinogenic amino acids
pyrrolysine:    PYL:O
selenocysteine: SEC:U
Glutamic acid or glutamine: GLX:Z
Unknown: XAA:X
Aspartic acid or Asparagine ASX:B
Fixed: 01/19/24
"""

# def main(argv):
def main():
    global hydrogen_B
    global lines_total
    global lines_total2
    global refined
    global HC_total2
    global VDW_total
    global repulsive
    global london
    global surface_tension
    global hydrophobicity
    global ligand_name
    ligand_name = "l_result.tmp" 

    # ****************************************************************************************************** #
    mapVolume = {}
    mapMass = {}
    mapHydro = {}

    seq = ''
    length = 0.0
    Nonpolar_Aliphatic = 0.0
    Aromatic = 0.0
    Polar_Uncharged = 0.0
    Positively_Charged = 0.0
    Negatively_Charged  = 0.0
    mass = 0.0
    volume = 0.0
    hydro = 0.0

    # ****************************************************************************************************** #
    global AA
    AA = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I","L", "K", "M", "F", "P", "S", "T", "W", "Y", "V","O","U","B","Z","X","J"]

    global AADipeptide
    AADipeptide = {}

    for i in AA:
        for j in AA:
            Dipeptide = i+j
            AADipeptide[Dipeptide] = int(0)

    global AATripeptide
    AATripeptide = {}

    for i in AA:
        for j in AA:
            for k in AA:
                kmer = (i+j+k)
                AATripeptide[kmer] = int(0)
    #if len(argv) < 4:
    #    sys.exit(1)
    table = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C','GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P','SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'PYL': 'O', 'SEC': 'U','GLX':'Z','XAA':'X','ASX':'B','MSE':'J','UNK':'>'}

    aa_table = {}

    for t in table:
        aa_table[table[t]] = t
    
    import os
    import sys
    #dir_ = './data/'
    dir_ = sys.argv[1]
    pdbfiles = list()
    for file in os.listdir(dir_):
        if file.endswith(".pdb"):
            pdbfiles.append(os.path.join(dir_, file))

    f = open('features.txt','r')
    key_ = list()
    for l in f:
        key_ = l.split('\t')

    finaloutput = open("1_Reults.interaction_terms.csv","w")
    finaloutput.write(f"pdb\t"+'\t'.join(key_))
    finaloutput.write('\n')
    i = 1
    for pdb in pdbfiles:
        # 1) processa entreadas
        # [JC] Open PDF file
        atoms,seqaa = open_protein(pdb,table)
        if '>' not in seqaa:
            print (f"{i}) Processing {pdb}")
            protein_ligand,protein_bond,protein_tripos = open_protein_ligand(pdb.replace('pdb','pdb.mol2'))
        
            #print (f"2.1 PP_Distance: [JC] Hydrophobic contacts (HCs) / van der Waals interactions (VDWs) / Repulsive interactions (RIs)")
            # Protein-protein distance
            file1,feat = PP_Distance(atoms,protein_ligand)
            raiz_proteina(file1,atoms)
            
            # "Utiliza a sessao bonds do mol2 e o L_DIST.tmp (PL_DIST.tmp) pra extrair as raizes de apenas atomos que estao listados no PL_DIST.tmp
            raiz_ligante(protein_bond,protein_tripos)
        
            #print (f"2 Calculos")
            #print (f"2.1 P_centro_geom: [JC] Hydrophobic contacts (HCs) / van der Waals interactions (VDWs) / Repulsive interactions (RIs)")
            P_centro_geom()
            
            #input eh o L_ROOT_0.tmp gerado pelo raiz_ligante, e vai gerar L_ROOT.tmp
            L_centro_geom()
            
            #print (f"input dos angulos: outputs dos P/L_centro_geom()")
            hydrogen_B = angulos()
            
            """
            # 3) Saida
            ---> remove ---> salva_ligante() #####
            """
            saida_PDB(atoms) 
            salva_proteina()
            
            #print (f"4) refinamento do HB")
            #result_score_calc
            if hydrogen_B != 0:
                limit()
            
            feat[2] = calcula_RT(protein_ligand)
            #print (feat[2])
            #print (f"Calculating scores")
            aux,key_ = result_score_calc(pdb,pdb.replace('pdb','pdb.mol2'),feat,seqaa,hydrogen_B,aa_table)
    
            #print (f"Saving results")
            finaloutput.write(f"{pdb}\t" + '\t'.join(map(str, aux)) + '\n')
    
            #print (f"Removing temp files")
            deleta_temp()
            
            #import time
            #time.sleep(0.1)
            i+=1
        
    finaloutput.close()
    
    """
    
    # system("pause")
    """
    
if __name__ == "__main__":
    #main(sys.argv[1:])
    main()


# In[ ]:




