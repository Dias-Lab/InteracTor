#include <iostream>
#include <fstream>
#include <cstdio>
#include <math.h>
#include <cstring>
#include <string>
#include <stdio.h>
#include <cstdlib>
#include <stdlib.h>
#include <unistd.h>
#include <limits>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <map>
#define TAM 82
/*
 * GLM-Score v0.15.2 changes
 * contact surface and total surface tension calculation added
 * parameter files hydrofobicity.param and tension.param were added
 * using maping functions to load paramater files (more efficient)
*/
int hydrogen_B = 0;
int lines_total = 0;
int lines_total2= 0;
int refined = 0; 
double HC_total2 = 0;
double VDW_total = 0;
int repulsive = 0;
double london = 0;
double RT = 0;
double surface_tension = 0;
double hydrophobicity = 0;
/*
Accepted distances for different types of bonds
C-C	1.507
C=C	1.336
C-N	1.449
C=N	1.273
Accepted angles
C-C-C	112.4
C-N-C	121.4
C-C-N	111.2
*/
/*
VDW radii - protein
1.65 	N, NE, NH1, NH2, ND2, NE2, ND1
1.40	O, OD1, OD2, OE1, OE2, OH, 
1.85	SG
1.50	NZ
RESIDUE ATOM ALA  
ATOM  N   1.65         
ATOM  O   1.40
RESIDUE ATOM ARG 
ATOM  N   1.65         
ATOM  O   1.40         
ATOM  NE  1.65          
ATOM  NH1 1.65          
ATOM  NH2 1.65  
RESIDUE ATOM ASP  
ATOM  N   1.65         
ATOM  O   1.40         
ATOM  OD1 1.40          
ATOM  OD2 1.40   
RESIDUE ATOM ASN  
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  OD1 1.40          
ATOM  ND2 1.65    
RESIDUE ATOM CYS                
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  SG  1.85       
RESIDUE ATOM GLU 
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  OE1 1.40          
ATOM  OE2 1.40      
RESIDUE ATOM GLN                
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  OE1 1.40          
ATOM  NE2 1.65      
RESIDUE ATOM GLY  
ATOM  N   1.65          
ATOM  O   1.40       
RESIDUE ATOM HIS  
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  ND1 1.65          
ATOM  NE2 1.65  
RESIDUE ATOM ILE                
ATOM  N   1.65          
ATOM  O   1.40   
RESIDUE ATOM LEU  
ATOM  N   1.65          
ATOM  O   1.40     
RESIDUE ATOM LYS  
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  NZ  1.50    
RESIDUE ATOM MET                
ATOM  N   1.65          
ATOM  O   1.40    
RESIDUE ATOM PHE              
ATOM  N   1.65          
ATOM  O   1.40       
RESIDUE ATOM PRO               
ATOM  N   1.65         
ATOM  O   1.40      
RESIDUE ATOM SER   
ATOM  N   1.65         
ATOM  O   1.40         
ATOM  OG  1.40      
RESIDUE ATOM THR  
ATOM  N   1.65         
ATOM  O   1.40         
ATOM  OG1 1.40
RESIDUE ATOM TRP 
ATOM  N   1.65         
ATOM  O   1.40         
ATOM  NE1 1.65     
RESIDUE ATOM TYR 
ATOM  N   1.65          
ATOM  O   1.40          
ATOM  OH  1.40     
RESIDUE ATOM VAL                
ATOM  N   1.65          
ATOM  O   1.40     
RESIDUE ATOM ASX               
ATOM  N   1.65          
ATOM  O   1.40     
RESIDUE ATOM GLX                
ATOM  N   1.65          
ATOM  O   1.40
*/
/*
VDW radii - ligand
C.3	1.94
C.2	1.90
C.1	1.90
C.ar	1.85
O.3	1.74
O.2	1.66
O.w	1.77
N.3	1.87
N.2	1.86
N.1	1.86
N.ar	1.86
N.am	1.83
N.pl3	1.86
S.3	2.09
S.2	2.01
S.o	2.01
S.o2	2.01
F	1.77
Cl	2.00
Br	2.22
I	2.42
P	2.03
*/
//---------------------HB-acceptors-and-donnors---------------------------------------
//-------------protein-------------------------------------------------------
//Side chains
/*
Doadores	aa				aceitadores			aa
N(H2)H		lis				N				his, try	
N(H)H		arg				OCO				glutamic, aspartic acid
NH		arg				CO				glu, asp
NH		trp, pro, his	OH					        tyr, thr, ser
N(H)H		glu, asp		S				        met (nao considerar)
COH		tyr, thr, ser	anel (tyr, trp)	nao considerar
SH		cys
CH (todas)	nao considerar
asp = asparigina
glu = glutamina
*/
// Main chain
// N (grupo amino)= doador e aceitador
// O (carbonila)= aceitador
//-------------------------------------------------------------------------------
//---------------------------Ligand--------------------------------------------
/*
Donnors: S N O
Acceptors: S N O
*/
//----------------------------Root---------------------------------------------
/*
//----------------PROTEINA-------------------------------------------------
Main chain:
D/A				RD/RA
N (amino)		CA
O (carbonila)	C (da carbonila)
Side chains:
			D/A		RD/RA
All aa		        N		CA e C (do aa anterior) ****************
		O (aceitador)	C
GLY			--		--
LYS			NZ		CE
ARG			NE		CD CZ **********************
			NH1		CZ		
			NH2		CZ
HIS			ND1		CG CE1	*******************	
			NE2		CD2 CE1 *********************
ASP			OD1 (aceitador) CG
			OD2		CG
GLU			OE1 (aceitador) CG
			OE2		CG
ASN			OD1 (aceitador) CG
			ND2		CG
SER			OG		CB
THR			OG1		CB
TYR			OH		CZ
TRP			N??		???
CYS			S		CB
GLN			OE1 (aceitador)	CD
			NE2		CD
TRP			NE1		CD1 CE2 ***************************
PRO			N (so aceita)	CA CD C (aa anterior)*********************
MET			SD (so aceita) CG CE ************************
ABDGEZ
*/

using namespace std;
void split(const string& s, char c, vector<string>& v) // pronto
{

	string::size_type i = 0; 
	string::size_type j = s.find(c); 
	string a;
	while (j != string::npos) {
		a = s.substr(i, j-i);
		if(a != "\0"){
			//cout << a << "@" << endl;
			v.push_back(a);
		}
		//v.push_back(s.substr(i, j-i)); 
		i = ++j; 
		j = s.find(c, j); 
		if (j == string::npos){
		a = s.substr(i, s.length( ));
			if(a != "\0" && a != "\n"){
				//cout << a << "#" << endl;
				v.push_back(a);
			}
		}
		//v.push_back(s.substr(i, s.length( ))); 
   } 
}

map<string, double> create_map_tables(string param_file_name){ // pronto
    map<string,double> param_map;
    ifstream param_file(param_file_name.c_str());
    string param_line;
    vector<string> param_v;
    while(getline(param_file, param_line)){
        split(param_line, ':', param_v);
        param_map[param_v[0]] = atof(param_v[1].c_str());
        //cout << param_v[0] << ": " << param_v[1] << endl;
        param_v.clear();
    }    
    return param_map;    
}
//
void surface_tension_hydrophobicity_calculator(map<string,double> hydro_map, map<string,double> tension_map, string infile_name) // pronto
{
    ifstream infile;
   	infile.open(infile_name.c_str());
    vector<string> v;
    string p_dist_line;
    vector<string> aa_numbers;
    int counter = 0;
    bool exists = 0;
    while (getline(infile, p_dist_line)) {
        v.clear();
        split(p_dist_line, ' ', v);
        //cout << p_dist_line << endl;
        //cout << "v1 " << v[1] << endl;
        exists = find(aa_numbers.begin(), aa_numbers.end(), v[4].c_str()) != aa_numbers.end();
        aa_numbers.push_back(v[4]);
        counter++;
        if(exists == 0){
            hydrophobicity += hydro_map[v[2]];
            //cout << "v2\t" << v[2] << "\tvalue\t" << hydro_map[v[2]] << "\taanum\t" << v[4] << endl;
            surface_tension += tension_map[v[2]];
        }
    }        
}
//
bool fncomp (char lhs, char rhs) {return lhs<rhs;} // pronto
	struct classcomp {
		bool operator() (const char& lhs, const char& rhs) const {return lhs<rhs;}
};
//
string ReplaceAll(string str, const string& from, const string& to){
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos){
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
    return str;
}

char protein_name[100];
char ligand_name[100];
char ligand_pdb[100];
char ligand_type[100];
//ofstream result_score;

double summation(double a, double b){ //pronto
	int x=0;
	double sum=0;
	for(x=a;x<=b;x++){
		sum+=x;
	}
	return sum;
}

int open_ligand() // pronto
{
	int achou_mol2 = 0;
	ifstream inlig;
	inlig.open(ligand_name);
	if (! inlig){
			cout << "Error opening input file...\n\nPlease, be sure that the *.PDB file you want to open is placed in the same\ndirectory of this program, that you wrote the complete and correct name and\n.pdb extention, and try again.\n" << endl;
			system("pause");
			exit(1);
	} 
	int pos = 0;
	ofstream outlig;
	outlig.open("ligand.tmp");
	//char line[82];		 
	string line;
	while(!inlig.eof()) {
		rd: line[0]='*';
		// [JC] The file was created by X-TOOL 
		getline(inlig, line);
		if(line == "@<TRIPOS>ATOM") {
			pos = 0;
			line[0]=' ';
			while (line[0]==' ' && !inlig.eof()) {
				getline(inlig, line);
				if (line[0]!='@'){
					// [JC] ... get 52 characters (example)
					// [JC] 1 N          15.0880   10.7980   23.5470 N.4       1 MOL         0.2328
					outlig << line.substr(0,52) << "\n";  
				}
			} 
			achou_mol2 = 1;
		}		   
	}
	outlig << endl; //remove this once I fix the bug that makes it necessary to have an extra new line when reading the ligand input.
	if (achou_mol2==0){cout << "Input error: " << ligand_name << " doesn't look to be a *.mol2 file. " << "please check input files and order of input data in comand line.\n"; system("pause"); exit(1);}
	inlig.close();
	outlig.close();
}
//-----------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
void distancia() // Pronto
{
	ofstream pdistmp;
	pdistmp.open("P_DIST.tmp");
	ofstream ldistmp;
	ldistmp.open("L_DIST.tmp");
	ifstream ligtmp;
	ligtmp.open("ligand.tmp");
	//-----------------------CALCULOS----------------------------
	int n = 0;  //conta posicao do caractere
	int i = 0;
	//-----------------------------------------DADOS DO LIGANTE----------------------
	string atomo;
	string nome_atomo;
	string x;
	string y;
	string z;
	char enter[2];
	int num_atomo = 0;
	int num_atomo_copia = 0;
	double num_x;
	double num_y;
	double num_z;
	//--------------------------------------------------------------------------------------------
	//distancia permitida
	double permitido = 0;
	//RADII 1 E RADII 2
	double r1 = 0;
	double r2 = 0;
	//------------------------------------------DADOS DA PROTEINA----------------------------------
	string p_atomo;
	string p_nome_atomo;
	string p_aminoacido;
	string p_num_aminoacido;
	string p_x;
	string p_y;
	string p_z;
	string chain;
	int num_p_atomo = 0;
	int num_p_atomo_copia = 0;
	double num_px;
	double num_py;
	double num_pz;
	int counter;
	ifstream protmp;
	protmp.open("protein.tmp"); //vai abrir o arquivo temporario da proteina para poder encontrar a raiz do doador ou aceitador e doador ou aceitador
	ofstream RT_found;
	RT_found.open("RT_found.tmp");
	ofstream HC_lig_contact;
	HC_lig_contact.open("HC_lig_contact.tmp");
	ofstream HC_lig_contact_lim;
	HC_lig_contact_lim.open("HC_lig_contact_lim.tmp");
	//--------------------------------------------------------------------------------------------
	string line = "*";
	string hibr;
	vector<string> v;
	int found_hibr = 0;
	int found_p = 0;
	while (!ligtmp.eof()){
		found_hibr = 0;
		getline(ligtmp, line);
		//cout << line << endl;
		v.clear();
		line = ReplaceAll(line, "-", " -");
		split(line, ' ', v); 
		if (line.compare(0,1,"") == 0 || line.compare(0,1,"\n") == 0  || line.compare(0,1,"\0") == 0){break;}
		atomo = v[0];
		num_atomo = atoi(v[0].c_str());
		//	 cout << "convertendo caracteres do ligante\n" << endl;
		//num_x = atof(x);
		nome_atomo = v[1];
		num_x = atof(v[2].c_str());
		x = v[2];
		num_y = atof(v[3].c_str());
		y = v[3];
		//cout << num_y;
		//cout << endl;
		num_z = atof(v[4].c_str());
		z = v[4];
		//cout << num_z;
		//cout << nome_atomo << endl;
		hibr = v[5];
		//cout << hibr << endl;
		//--------LIGANTE-------------RADII-Raio-de-VDW-PERMITIDOS
		if(hibr.compare(0,3,"N.3") == 0){
			r2 = 1.87;found_hibr = 1;
			//cout << hibr << " " << r2 << endl; 
		}
		if(hibr.compare(0,3,"N.1") == 0 || hibr.compare(0,3,"N.2") == 0 || hibr.compare(0,4,"N.ar") == 0 || hibr.compare(0,5,"N.pl3") == 0){
			r2 = 1.86;found_hibr = 1;
			//cout << hibr << " " << r2 << endl; 
		}
		if(hibr.compare(0,3,"N.am") == 0){
			r2 = 1.83;found_hibr = 1;
			//cout << hibr << " " << r2 << endl; 
		}
		if(hibr.compare(0,3,"C.3") == 0){
			r2 = 1.94;found_hibr = 1;
			//cout << hibr << " " << r2 << endl; 
		}
		if(hibr.compare(0,3,"C.1") == 0 || hibr.compare(0,3,"C.2") == 0){
			r2 = 1.90;found_hibr = 1;
		}
		if(hibr.compare(0,4,"C.ar") == 0){
			r2 = 1.85;found_hibr = 1;
		}
		if(hibr.compare(0,3,"O.3") == 0){
			r2 = 1.74;found_hibr = 1;
		}
		if(hibr.compare(0,3,"O.2") == 0 || hibr.compare(0,5,"O.co2") == 0){ //check the radius of C.co2
			r2 = 1.66;found_hibr = 1;
			//cout << "hibr: " << hibr << " " << r2 << endl; 
		}
		if(hibr.compare(0,1,"F") == 0 || hibr.compare(0,3,"O.w") == 0){
			r2 = 1.77;found_hibr = 1;
		}
		if(hibr.compare(0,1,"S") == 0){
			r2 = 2.01;found_hibr = 1;
		}
		if(hibr.compare(0,3,"S.3") == 0){
			r2 = 2.09;found_hibr = 1;
		}
		if(hibr.compare(0,2,"Cl") == 0){
			r2 = 2.00;found_hibr = 1;
		}
		if(hibr.compare(0,2,"Br") == 0){
			r2 = 2.22;found_hibr = 1;
		}
		if(hibr.compare(0,1,"I") == 0){
			r2 = 2.42;found_hibr = 1;
		}
		if(hibr.compare(0,1,"P") == 0){
			r2 = 2.03;found_hibr = 1;
		}
		//--------LIGANTE-------------RADII-Raio-de-VDW-PERMITIDOS
		//cout << "*" << atomo << "*" << espaco_1 << "*" << nome_atomo << "*" << espaco_2 << "*" << espaco_3 << "*" << x << "*" << y << "*" << z << "*" << espaco_5 << "*" << hibr << "*" << enter;
		//cout << endl;
		//cout << "number: " << num_x << "name: " << v[2] << endl;
		double DIS2 = 0;
		string linep;
		while (!protmp.eof()){
			found_p = 0;
			getline(protmp, linep);
			//cout << linep << endl;
			v.clear();
			line = ReplaceAll(line, "-", " -");
			split(linep, ' ', v);
			if(v[0] == "" || v[0] == "\0" || v[0] == "\n"){break;}
			//cout << "*" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << px << "*" << py << "*" << pz << "*" << p_enter;
			if (linep.compare(0,1,"") == 0 || linep.compare(0,1,"\n") == 0 || linep.compare(0,1,"\0") == 0){break;}
			//cout << endl;
			//cout << "convertendo caracteres da proteina\n" << endl;
			p_atomo = v[0];
			num_p_atomo = atoi(v[0].c_str());
			p_nome_atomo = v[1];
			p_aminoacido = v[2];
			chain = v[3];
			p_num_aminoacido = v[4];
			num_px = atof(v[5].c_str());
			p_x = v[5];
			//cout << num_px << endl;
			num_py = atof(v[6].c_str());
			p_y = v[6];
			//cout << num_py << endl;
			num_pz = atof(v[7].c_str());
			p_z = v[7];
			//cout << num_pz << endl;
			//double DDA = sqrt(((num_x-num_px)*(num_x-num_px))+((num_y-num_py)*(num_y-num_py))+((num_z-num_pz)+(num_z-num_pz))); 
			//distancia entre o atomo do doador e aceitador  
			//angulos
			//cos-1=[(AB)2 + (AC)2 - (BC)2 / 2(AB)*(AC)]
			double DX = (num_x - num_px);
			double DY = (num_y - num_py);
			double DZ = (num_z - num_pz);
			double EDX = (DX*DX);
			double EDY = (DY*DY);
			double EDZ = (DZ*DZ);
			double SD = (EDX+EDY+EDZ);
			double DIS = sqrt(SD);
			//cout << DIS;
			DIS2 = DIS;
			//cout << endl;
			//CONTINUACAO ACHA DISTANCIA PERMITIDA
			//cout << endl << "*" << p_nome_atomo << "*" << endl;
			//cout << endl << "*" << nome_atomo << "*" << endl;
			if (p_nome_atomo.compare(0,1,"N") == 0){
				r1 = 1.65; found_p = 1;
				//cout << " found : " << p_nome_atomo << endl;
			}
			if (p_nome_atomo.compare(0,1,"O") == 0){
				r1 = 1.40;found_p = 1;
			}
			if (p_nome_atomo.compare(0,2,"NZ") == 0){
				r1 = 1.50;found_p = 1;
			}
			if (p_nome_atomo.compare(0,1,"S") == 0){
				r1 = 1.85;found_p = 1;
			}
			if (p_nome_atomo.compare(0,1,"C") == 0){
				r1 = 1.76; found_p = 2;
			}
			if (p_nome_atomo.compare(0,2,"CA") == 0 || p_nome_atomo.compare(0,2,"CB") == 0 || p_nome_atomo.compare(0,2,"CD") == 0 || p_nome_atomo.compare(0,2,"CE") == 0 || p_nome_atomo.compare(0,2,"CG") == 0 || p_nome_atomo.compare(0,2,"CH") == 0 || p_nome_atomo.compare(0,2,"CZ") == 0){
				r1 = 1.87; found_p = 2;
			}
			/*
			1.65 	N, NE, NH1, NH2, ND2, NE2, ND1
			1.40	O, OD1, OD2, OE1, OE2, OH, 
			1.85	SG
			1.50	NZ
			*/
			permitido = r1+r2;
			double permitido1 = permitido-0.7; 
			double extended = 0.7;
			// RT
			if (DIS2<=(permitido + extended) && nome_atomo.compare(0,1,"H") != 0  && found_hibr == 1 && found_p == 1){
				//cout << DIS2 << " " << p_nome_atomo << " " << atomo << " " << hibr  << " r1 " << r1 << " r2 " << r2 << endl; 
				//cout << atomo[1] << atomo[2] << atomo[3] << atomo[4] << atomo[5] << atomo[6] << endl;
				for(i=0; i<(6 - atomo.size()); i++){RT_found << " "; }
				RT_found << atomo << endl;
			}
			// VDW
			// [JC] van der Waals interactions (VDWs)
			if (hibr.compare(0,1,"H") != 0 && atomo.compare(0,1,"H") != 0  && found_hibr == 1 && found_p != 0){
				//funcao para o VDW
				double VDW_radii = r1+r2;
				//cout << VDW_radii << endl;
				double var_0 = VDW_radii/DIS2;
				//cout << r1 << " " << r2 << endl;
				//cout << var_0 << endl;
				double exp1 = 8.0;
				double exp2 = 4.0;
				double var_1 = pow(var_0, exp1); 
				double var_2 = pow(var_0, exp2);
				double VDW_int = var_1-(2*(var_2));
				if (VDW_int>=100){VDW_int = 100;}
				if (VDW_int<=100){VDW_total += VDW_int;} 
				VDW_int = 0;
				//cout << "VDW: " << VDW_total << "r1: " << r1 << "r2: "<< r2 << endl;
			}//else{cout << "excluded: " << hibr << ", " << nome_atomo << endl;}//if hibr 0 != H
			// HC
			// [JC] Hydrophobic contacts (HCs) 
			double HC = 0;
			double HC2 = 0;
			double HC_exp = 0;
			if (p_nome_atomo.compare(0,1,"C") == 0  && nome_atomo.compare(0,1,"C") == 0 && found_hibr == 1 && found_p != 0) {
				double HC_VDW = permitido;
				double HC_allowed_1 = (HC_VDW+0.5);
				double HC_allowed_2 = (HC_VDW+2.0);
				//---------------------HC2--------------------------------
				if(DIS2<=HC_allowed_1){HC2 = 1;}
				if((DIS2>HC_allowed_1) && (DIS2<=HC_allowed_2)){
					HC2 = ((1/1.5)*(pow((HC_VDW+2.0),2)-pow(DIS2,2)));
					//cout << "HC_total: " << HC_total2 << " HC2: " << HC2 << " perimitido1: " << permitido1 << " r1: " << r1 << " r2: " << r2 << endl;
					//cout << "*" << HC << endl;
				}
				if(DIS2>HC_allowed_2){HC2 = 0;}
				//if(DIS2>HC_allowed_2){HC_total2+=0;}
				HC_total2+=HC2;
				//cout << "HC_total: " << HC_total2 << " HC2: " << HC2 << " perimitido1: " << permitido1 << " r1: " << r1 << " r2: " << r2 << endl;
			}
			//----------------------- REPULSIVE INTERACTIONS DIS2 <= permitido
			// [JC] Repulsive interactions (RIs)
			if (DIS2<=permitido && nome_atomo.compare(0,1,"H") != 0  && found_hibr == 1) {
				repulsive++;
				if (london+1/DIS2 != numeric_limits<double>::infinity( )){
					london+=1/DIS2;
				}
			}
			//---------------------HC2--------------------------------
			// [JC] Hydrophobic contacts (HCs) 2
			if (/*DIS2 <= permitido1 ||*/DIS2>=permitido1 && DIS2<=permitido  && found_hibr == 1 && found_p == 1){
				if ((p_nome_atomo.compare(0,1,"N") == 0 || p_nome_atomo.compare(0,1,"O") == 0 || p_nome_atomo.compare(0,1,"S") == 0) && (nome_atomo.compare(0,1,"N") == 0 || nome_atomo.compare(0,1,"O") == 0 || nome_atomo.compare(0,1,"S") == 0)){
					//cout << "DISTANCIA: " << DIS2 << endl;
					//cout << p_nome_atomo << "  r1: " << r1 << endl;
					//cout << hibr << " r2: " << r2 << endl;
					//cout << "permitido (r1+r2): " << permitido << endl;
					//cout << "(permitido - 0.7): " << permitido1 << endl << endl;
					//cout << "r1+r2 = " << r1+r2 << endl;
					n = 0;
					for(n=0; n<(7-atomo.size()); n++){ldistmp << " ";}
					ldistmp << atomo << "  " << nome_atomo;
					for(n=0; n<(7-nome_atomo.size()); n++){ldistmp << " ";}
					for(n=0; n<(10-x.size()); n++){ldistmp << " ";}
					ldistmp << x;
					for(n=0; n<(10-y.size()); n++){ldistmp << " ";}
					ldistmp << y;
					for(n=0;n<(10-z.size()); n++){ldistmp << " ";}
					ldistmp << z << " " << hibr;
					for(n=0; n<(5-hibr.size()); n++){ldistmp << " ";}
					ldistmp << endl;
					n = 0;
					//proteina: salva arquivo tmp com os atomos mais pr�imos
					if (p_atomo.size() < 6){
						for(n=0; n<(6-p_atomo.size()); n++){pdistmp << " ";}
					}
					pdistmp << p_atomo << " " << p_nome_atomo;
					if (p_nome_atomo.size() < 4){
						for(n=0; n<(4-p_nome_atomo.size()); n++){pdistmp << " ";}
					}
					pdistmp << p_aminoacido;
					if (p_aminoacido.size() < 4){
						for(n=0; n<(4-p_aminoacido.size()); n++){pdistmp << " ";}
					}
					pdistmp << chain;
					if (p_num_aminoacido.size() < 4){
						for(n=0; n<(4-p_num_aminoacido.size()); n++){pdistmp << " "; }
					}
					if (p_num_aminoacido.size() >= 4 && p_num_aminoacido.size() < 9){
						for(n=0; n<(9-p_num_aminoacido.size()); n++){pdistmp << " "; }
					}
					pdistmp << p_num_aminoacido;
					if (p_x.size() < 9){
						for(n=0; n<(9-p_x.size()); n++){pdistmp << " ";}
					}
					pdistmp << p_x;
					if (p_y.size() <= 8){
						for(n=0; n<(8-p_y.size()); n++){pdistmp << " ";}
					}
					pdistmp << p_y;
					if (p_z.size() <= 8){
						for(n=0; n<(8-p_z.size()); n++){pdistmp << " ";}
						pdistmp << p_z << " \n";
					}else{pdistmp << p_z << " \n";}
					n = 0;
					//cout << "distancia " << DIS2 << endl;
				}
				//}
			}//if (/*DIS2 <= permitido1 ||*/DIS2>=permitido1 && DIS2<=permitido  && found_hibr == 1){
		}
		protmp.clear();              // forget we hit the end of file
		protmp.seekg(0, ios::beg);   // move to the start of the file
	}
	protmp.close(); 
	pdistmp.close();
	ldistmp.close();  
	RT_found.close();
	HC_lig_contact.close();
	HC_lig_contact_lim.close();
}
//-----------------------------------------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------
int open_protein() // Pronto
{
	int achou_pdb = 0;
	ifstream inprot;
	inprot.open(protein_name);
	if (! inprot)
	{
			cout << "Error opening input file...\n\nPlease, be sure that the *.PDB file you want to open is placed in the same\ndirectory of this program, that you wrote the complete and correct name and\n.pdb extention, and try again.\n" << endl;
			exit(1);
	} 
	ofstream outprot;
	outprot.open("protein.tmp");
	// [RS] linha que o programa le do arquivo PDB (Example)
	// [JC] ATOM      1  N   PRO A   2      31.242   3.064  39.284  1.00 39.90           N
	char line[82];		 
	int lx = 0;
	while (lx!=80)
	{
		line[lx]==' '; lx++;
	}
	//lx = 4;
	lx = 0;
	int lx1 = 1;
	char ch[2];
	while(!inprot.eof())
	{
		/*
		if (line[61]=='\n')
		{inprot.ignore(62, '\n');}
		if (line[62]=='\n')
		{inprot.ignore(63, '\n');}
		if (line[63]=='\n')
		{inprot.ignore(64, '\n');}
		if (line[64]=='\n')
		{inprot.ignore(65, '\n');}
		if (line[65]=='\n')
		{inprot.ignore(66, '\n');}
		//------------------------------
		*/
		// [JC] Open the PDB file
		inprot.get(ch[0]);
		line[0]=ch[0];
		//cout << "*" << ch[0] << "*" << endl;
		while (!inprot.eof() && line[lx]!='\n')
		{
			inprot.get(ch[0]);
			line[lx1]=ch[0];
			lx++; lx1++;
			//cout << "*" << ch[0] << line[lx] << "*" << endl;
		}
		//inprot.get(line,sizeof(line),'\0');
		line[80]='\n';
		lx = 0;
		lx1 = 1;
		//cout << "line 0 =" << line[0] << "*" << endl;
		//cout << "line80 =" << line[80] << "*" << endl << endl;
		//if (line[80]!='\n')
		//{inprot.ignore (1, '\n');}
		//{inprot.get(line[70]);}
		// [JC] get ATOM line
		// [JC] ATOM      1  N   PRO A   2      31.242   3.064  39.284  1.00 39.90           N
		if (line[0]=='A' && line[1]=='T' && line[2]=='O' && line[3]=='M' && line[80]=='\n')
		{	
			// Why ?   
			lx = 4;
			while (lx!=5) 
			{lx++;} 
			// Why ?
			while (lx!=11)
			{outprot.put(line[lx]); lx++;}

			// [JC] get Element (Third column)
			outprot.put(' '); 
			outprot.put(line[13]);
			outprot.put(line[14]);
			outprot.put(line[15]);
			// [JC] get chain/ strand ID (fourth column)
			outprot.put(' ');
			outprot.put(line[17]);
			outprot.put(line[18]);
			outprot.put(line[19]);
			// [JC] get residuo number
			outprot.put(' ');
			outprot.put(line[21]);		
			outprot.put(line[22]);
			outprot.put(line[23]);
			outprot.put(line[24]);
			outprot.put(line[25]);
			outprot.put(' ');

			while (lx!=30){lx++;}
			while (lx!=54)
			{
				outprot.put(line[lx]); lx++;
			}
			outprot.put(' ');
			while (lx!=0)
			{lx--;}
			outprot.put('\n');
			achou_pdb = 1;
		}
		while (lx!=80)
		{line[lx]==' '; lx++;}
		lx = 0;
	}
	if (achou_pdb==0){cout << "Input error: " << protein_name << " doesn't look to be a *.pdb file. " << "please check input files and order of input data in comand line.\n"; system("pause"); exit(1);}
	inprot.close();
	outprot.close();
}
//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
void raiz_proteina()
{
	ofstream root;
	root.open("P_ROOT_0.tmp");
	ifstream AD;
	char p_atomo[8];
	char p_nome_atomo[5];
	char p_aminoacido[4];
	char p_num_aminoacido[8];
	char px[9];
	char py[9];
	char pz[9];
	char p_enter[3];
	int num_p_num_aminoacido;
	int num_pp_num_aminoacido;
	int num_p_atomo = 0;
	int num_p_atomo_copia = 0;
	double num_px;
	double num_py;
	double num_pz;
	int n = 0;
	//copia dados da proteina 2 P_DIST (atomos com a distancia aceita)
	char pp_atomo[8];
	char pp_nome_atomo[5];
	char pp_aminoacido[4];
	char pp_num_aminoacido[8];
	char ppx[9];
	char ppy[9];
	char ppz[9];
	char pp_enter[3];
	int num_pp_atomo = 0;
	int num_pp_atomo_copia = 0;
	double num_ppx;
	double num_ppy;
	double num_ppz;
	AD.open("P_DIST.tmp");
	// tentar abrir arquivo da proteina
	while (!AD.eof())
	{
		AD.get(pp_atomo,sizeof(pp_atomo),'\0');  
		AD.get(pp_nome_atomo,sizeof(pp_nome_atomo),'\0');  
		AD.get(pp_aminoacido,sizeof(pp_aminoacido),'\0'); 
		AD.get(pp_num_aminoacido,sizeof(pp_num_aminoacido),'\0');  
		AD.get(ppx,sizeof(ppx),'\0');  
		AD.get(ppy,sizeof(ppy),'\0');  
		AD.get(ppz,sizeof(ppz),'\0');  
		AD.get(pp_enter,sizeof(pp_enter),'\0');  
		//cout << "1" << pp_atomo <<  "2" << pp_nome_atomo <<  "3" << pp_aminoacido <<  "4" << pp_num_aminoacido <<  "5" << ppx <<  "6" << ppy <<  "7" << ppz <<  "8" << pp_enter << endl;
		//cout << pp_nome_atomo[0] << pp_nome_atomo[1] << pp_nome_atomo[2] << endl;
		ifstream AD_root;
		AD_root.open("protein.tmp");
		int aa = 0;
		while (!AD_root.eof())
		{
			aa++;
			AD_root.get(p_atomo,sizeof(p_atomo),'\0');  
			AD_root.get(p_nome_atomo,sizeof(p_nome_atomo),'\0');  
			AD_root.get(p_aminoacido,sizeof(p_aminoacido),'\0'); 
			AD_root.get(p_num_aminoacido,sizeof(p_num_aminoacido),'\0');  
			AD_root.get(px,sizeof(px),'\0');  
			AD_root.get(py,sizeof(py),'\0');  
			AD_root.get(pz,sizeof(pz),'\0');  
			AD_root.get(p_enter,sizeof(p_enter),'\0');  
			// cout << "1" << p_atomo <<  "2" << p_nome_atomo <<  "3" << p_aminoacido <<  "4" << p_num_aminoacido <<  "5" << px <<  "6" << py <<  "7" << pz <<  "8" << p_enter << endl;
			// -------------------------------N-PROLINA
			// [JC] 
			if (pp_nome_atomo[0]=='N' && pp_aminoacido[0]=='P' && pp_aminoacido[1]=='R' && pp_aminoacido[2]=='O')
			{
				if (aa!=1)
				{	
					num_p_num_aminoacido = atoi(p_num_aminoacido);
					num_pp_num_aminoacido = atoi(pp_num_aminoacido);
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]==' ' && num_p_num_aminoacido==(num_pp_num_aminoacido -1 ))
					{
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;  			
						//------cria arquivo temp com as raizes
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;  			
					}//num_p_num_aminoacido==(num_pp_num_aminoacido-1))
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='A' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}		
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}		
				}//if (aa!=1)
				if (aa==1) //se nao for o primeiro aminoacido
				{	
					// Esse if estava sem
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='A' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7]){							
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
						//cout << "1" << p_atomo <<  "2" << p_nome_atomo <<  "3" << p_aminoacido <<  "4" << p_num_aminoacido <<  "5" << px <<  "6" << py <<  "7" << pz <<  "8" << p_enter << endl;  
						//cout << "1" << pp_atomo <<  "2" << pp_nome_atomo <<  "3" << pp_aminoacido <<  "4" << pp_num_aminoacido <<  "5" << ppx <<  "6" << ppy <<  "7" << ppz <<  "8" << pp_enter << endl;	   
					//if ........
					}	
				}//if (aa==1)
			}
			//----------------------------------N-TODAS-PROTEINAS-EXCETO-PRO
			if (pp_nome_atomo[0]=='N' && (pp_aminoacido[0]!='P' || pp_aminoacido[1]!='R' || pp_aminoacido[2]!='O'))
			{
				//----------------------------------N
				if (pp_nome_atomo[1]==' ')
				{
					if (aa!=1)//se nao for o prmeiro aminoacido
					{	
						num_p_num_aminoacido = atoi(p_num_aminoacido);
						num_pp_num_aminoacido = atoi(pp_num_aminoacido);
						int num_pp_num_aminoacido_2 =(num_pp_num_aminoacido-1);
						if (p_nome_atomo[0]=='C' && p_nome_atomo[1]==' ' && num_p_num_aminoacido==num_pp_num_aminoacido_2)
						{
							//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter << endl;  			
							root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
						}//num_p_num_aminoacido==(num_pp_num_aminoacido -1 ))
						if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='A' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
						{	
							//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
							root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
						}		
					}//if (aa!=1)
					if (aa==1)//se for o primeiro amino�ido
						{	
						if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='A' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
						{	
							//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
							root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
							//cout << "1" << pp_atomo <<  "2" << pp_nome_atomo <<  "3" << pp_aminoacido <<  "4" << pp_num_aminoacido <<  "5" << ppx <<  "6" << ppy <<  "7" << ppz <<  "8" << pp_enter << endl;	   
						}//if ........
					}//if (aa==1)
				}//if (pp_nome_atomo[1]==' ')
				
				//----------------------------------NZ-LYS
				if (pp_nome_atomo[1]=='Z' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='L' && pp_aminoacido[1]=='Y' && pp_aminoacido[2]=='S')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='E' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='Z')
				
				//----------------------------------NE-ARG
				if (pp_nome_atomo[1]=='E' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='A' && pp_aminoacido[1]=='R' && pp_aminoacido[2]=='G')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='Z' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='E')
				
				//----------------------------------NH1/2-ARG
				if (pp_nome_atomo[1]=='H' && (pp_nome_atomo[2]=='1' || pp_nome_atomo[2]=='2') && pp_aminoacido[0]=='A' && pp_aminoacido[1]=='R' && pp_aminoacido[2]=='G')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='Z' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						//------cria arquivo temp com as raizes
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='H')
				
				//----------------------------------ND1-HIS
				if (pp_nome_atomo[1]=='D' && pp_nome_atomo[2]=='1' && pp_aminoacido[0]=='H' && pp_aminoacido[1]=='I' && pp_aminoacido[2]=='S')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='E' && p_nome_atomo[2]=='1' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])   
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='D')
				
				//----------------------------------NE2-HIS
				if (pp_nome_atomo[1]=='E' && pp_nome_atomo[2]=='2' && pp_aminoacido[0]=='H' && pp_aminoacido[1]=='I' && pp_aminoacido[2]=='S'){
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[2]=='2' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='E' && p_nome_atomo[2]=='1' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='E')
				
				//---------------------------ND2--ASN	
				if (pp_nome_atomo[1]=='D' && pp_nome_atomo[2]=='2' && pp_aminoacido[0]=='A' && pp_aminoacido[1]=='S' && pp_aminoacido[2]=='N')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='D')
				//---------------------------NE1--GLN	
				if (pp_nome_atomo[1]=='E' && pp_nome_atomo[2]=='2' && pp_aminoacido[0]=='G' && pp_aminoacido[1]=='L' && pp_aminoacido[2]=='N')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='E')
				
				//---------------------------NE1--TRP	
				if (pp_nome_atomo[1]=='E' && pp_nome_atomo[2]=='1' && pp_aminoacido[0]=='T' && pp_aminoacido[1]=='R' && pp_aminoacido[2]=='P')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[2]=='1' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='E' && p_nome_atomo[3]=='2' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........
				}//if (pp_nome_atomo[1]=='E')
			}//p_nome_atomo[0]=='N'
			// Ate aqui ()
			if (pp_nome_atomo[0]=='O')
			{
				//----------------------------------O
				if (pp_nome_atomo[1]==' ')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]==' ' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]==' ')
				//---------------------------OD1/2--ASP	
				if (pp_nome_atomo[1]=='D' && (pp_nome_atomo[2]=='1' || pp_nome_atomo[2]=='2') && pp_aminoacido[0]=='A' && pp_aminoacido[1]=='S' && pp_aminoacido[2]=='P')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='D')
				//---------------------------OE1/2--GLU	
				if (pp_nome_atomo[1]=='E' && (pp_nome_atomo[2]=='1' || pp_nome_atomo[2]=='2') && pp_aminoacido[0]=='G' && pp_aminoacido[1]=='L' && pp_aminoacido[2]=='U')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='E')
				//---------------------------OD1--ASN	
				if (pp_nome_atomo[1]=='D' && pp_nome_atomo[2]=='1' && pp_aminoacido[0]=='A' && pp_aminoacido[1]=='S' && pp_aminoacido[2]=='N')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[2]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7]){
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='D')
				//---------------------------OG--SER	
				if (pp_nome_atomo[1]=='G' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='S' && pp_aminoacido[1]=='E' && pp_aminoacido[2]=='R')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='B' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
					//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
					root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='G')
				// AQUI
				//---------------------------OG1--THR	
				if (pp_nome_atomo[1]=='G' && pp_nome_atomo[2]=='1' && pp_aminoacido[0]=='T' && pp_aminoacido[1]=='H' && pp_aminoacido[2]=='R')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='B' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7]){	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='G')
				//---------------------------OH--TYR	
				if (pp_nome_atomo[1]=='H' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='T' && pp_aminoacido[1]=='Y' && pp_aminoacido[2]=='R')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='Z' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='G')
				

				//---------------------------OE1--GLN	
				if (pp_nome_atomo[1]=='E' && pp_nome_atomo[2]=='1' && pp_aminoacido[0]=='G' && pp_aminoacido[1]=='L' && pp_aminoacido[2]=='N')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='D' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
					{	
						//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
						root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					}//if ........	
				}//if (pp_nome_atomo[1]=='E')

				//---------------------------OXT--TODOS		
				if (pp_nome_atomo[1]=='X' && pp_nome_atomo[2]=='T')
				{
					if (p_nome_atomo[0]=='C' && p_nome_atomo[1]==' ' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
						{	
							//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
							root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
						}//if ........	
				}//if (pp_nome_atomo[1]=='X' && pp_nome_atomo[2]=='T')
			}//p_nome_atomo[0]=='O'

			//-------------------------------------S
			//--------------------------------SD-MET
			if (pp_nome_atomo[0]=='S' && pp_nome_atomo[1]=='D' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='M' && pp_aminoacido[1]=='E' && pp_aminoacido[2]=='T')
			{
				if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='E' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
				{	
				//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
				root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
				}//if ........		  
				if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='G' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
				{	
					//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
					root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
				}//if ........		  
			}//if (pp_nome_atomo[0]=='S')

			if (pp_nome_atomo[0]=='S' && (pp_nome_atomo[1]==' ' || pp_nome_atomo[1]=='G') && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='C' && pp_aminoacido[1]=='Y' && pp_aminoacido[2]=='S')
			{
				if (p_nome_atomo[0]=='C' && p_nome_atomo[1]=='B' && p_nome_atomo[3]==' ' && pp_aminoacido[0]==p_aminoacido[0] && pp_aminoacido[1]==p_aminoacido[1] && pp_aminoacido[2]==p_aminoacido[2] && pp_num_aminoacido[0]==p_num_aminoacido[0] && pp_num_aminoacido[1]==p_num_aminoacido[1] && pp_num_aminoacido[2]==p_num_aminoacido[2] && pp_num_aminoacido[3]==p_num_aminoacido[3] && pp_num_aminoacido[4]==p_num_aminoacido[4] && pp_num_aminoacido[5]==p_num_aminoacido[5] && pp_num_aminoacido[6]==p_num_aminoacido[6] && pp_num_aminoacido[7]==p_num_aminoacido[7])
				{	
					//cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;   
					root << "R" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
				}//if ........		  
			}//if (pp_nome_atomo[0]=='S')
		}//while(!AD_root.eof())
	//cout << "B" << pp_atomo <<  "*" << pp_nome_atomo <<  "*" << pp_aminoacido <<  "*" << pp_num_aminoacido <<  "*" << ppx <<  "*" << ppy <<  "*" << ppz <<  "*" << pp_enter << endl;	   
	root << "B" << pp_atomo << pp_nome_atomo << pp_aminoacido << pp_num_aminoacido << ppx << ppy << ppz << pp_enter;
	AD_root.close();
	}//while (!AD.eof())
	AD.close();
	root.close();

}

/////////////////////////////////////////////////////////////////////////////////////////////////
//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
void raiz_ligante()
{

	char line[82];		 //linha que o programa le do arquivo PDB
	char ordem[30];
	char atomo_1[30];
	char atomo_2[30];
	int num_ordem = 0;
	int num_atomo_1 = 0;
	int num_atomo_2 = 0; 
	ofstream bond_type;
	bond_type.open("limit_type.tmp");
	ofstream bonds;
	bonds.open("bonds.tmp");
	ifstream LAD; //ARQUIVO MOL2: ACHA ATOMO LIGANTE (ACEITADOR / DOADOR)
	LAD.open(ligand_name);
	while (!LAD.eof())
	{
	
		rd: line[0]='*';
		while (line[0]!='@' && !LAD.eof()){LAD.get(line[0]);}
		int pos = 1; 
		if (line[0]=='@')
		{   
			while(pos!=14) 
			{LAD.get(line[pos]); pos++;}
		}//if (line[0]=='@')
		if (line[0]=='@' && line[1]=='<' && line[2]=='T' && line[3]=='R' && line[4]=='I' && line[5]=='P' && line[6]=='O' && line[7]=='S' && line[8]=='>' && line[9]!='B' && line[10]!='O' && line[11]!='N' && line[12]!='D' && line[13]!='\n')
		{
			//pos = 0;
			goto rd;	
		}//!BOND
		if (line[0]=='@' && line[1]=='<' && line[2]=='T' && line[3]=='R' && line[4]=='I' && line[5]=='P' && line[6]=='O' && line[7]=='S' && line[8]=='>' && line[9]=='B' && line[10]=='O' && line[11]=='N' && line[12]=='D' && line[13]=='\n')
		{    
			// pos = 0;
			line[0]=' ';
			//contador de coluna
			point:
			while (ordem[0]!='@' && !LAD.eof())
			{
				int c1 = 0;
				int c2 = 0;
				ordem[0] = ' ';
				while (ordem[0]==' '  && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(ordem[0]);
				}
				c1 = 0;
				c2 = 1;
				while (ordem[c1]!=' ' && ordem[c1]!='@'  && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(ordem[c2]);
					c1++;
					c2++;
				}
				//cout << ordem << endl;
				atomo_1[0] = ' ';
				while (atomo_1[0]==' '  && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(atomo_1[0]);
				}
				c1 = 0;
				c2 = 1;
				while (atomo_1[c1]!=' ' && atomo_1[c1]!='@' && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(atomo_1[c2]); //cout << atomo_1[c1];
					c1++;
					c2++;
				}
				atomo_2[0] = ' ';
				while (atomo_2[0]==' ' && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(atomo_2[0]);
				}
				c1 = 0;
				c2 = 1;
				while (atomo_2[c1]!=' ' && atomo_2[c1]!='@'  && !LAD.eof())
				{//LAD.get(ordem, sizeof(ordem),'\0');
					LAD.get(atomo_2[c2]); //cout << atomo_2[c1];
					c1++;
					c2++;
				}
				char tipo [8];
				tipo[0] = ' ';
				tipo[1] = ' ';
				tipo[2] = ' ';
				tipo[3] = ' ';
				tipo[4] = ' ';
				tipo[5] = ' ';
				tipo[6] = ' ';
				tipo[7] = ' ';
				if (atomo_2[c1]==' ' && atomo_2[c1]!='@'  && !LAD.eof())
				{
					while (tipo[0]==' '){LAD.get(tipo[0]);}
				}
				if (tipo[0]!=' '){
					int ct1 = 0; 
					int ct2 = 1; 
					while(tipo[ct1]!='\n'){LAD.get(tipo[ct2]); ct1++; ct2++;}
				}
				//cout << tipo << endl;
				if (ordem[0]!='@' && ordem[c1]!='@' && atomo_1[c1]!='@' && atomo_2[c1]!='@' && atomo_1[c2]!='@' && atomo_2[c2]!='@')
				{
					///////
					num_ordem = atoi(ordem);
					num_atomo_1 = atoi(atomo_1); 
					num_atomo_2 = atoi(atomo_2); 
					//cout <<  "ordem*" << ordem  << "*ordem" << atomo_1 <<  "atomo1*" << atomo_2 <<  "atomo2*" << endl;
					//cout << "numeros:" << endl;
					// cout << "*" << num_ordem << "*" << " " << "*" << num_atomo_1 << "*" << " " << "*" << num_atomo_2 << "*" << endl;
					// bonds << num_atomo_1 << num_atomo_2 << endl;
					if (num_atomo_1<=9 && num_atomo_1!=0)
					{bonds << "     " << num_atomo_1; bond_type << "     " << num_atomo_1 << endl;}
					if (num_atomo_1<=99 && num_atomo_1>=10)
					{bonds << "    " << num_atomo_1; bond_type << "    " << num_atomo_1 << endl;}
					if (num_atomo_1<=999 && num_atomo_1>=100)
					{bonds << "   " << num_atomo_1; bond_type << "   " << num_atomo_1 << endl;}
					if (num_atomo_1<=9999 && num_atomo_1>=1000)
					{bonds << "  " << num_atomo_1; bond_type << "  " << num_atomo_1 << endl;}
					if (num_atomo_1<=99999 && num_atomo_1>=10000)
					{bonds << " " << num_atomo_1; bond_type << " " << num_atomo_1 << endl;}
					if (num_atomo_1<=999999 && num_atomo_1>=100000)
					{bonds << num_atomo_1; bond_type << num_atomo_1 << endl;}
					//cout << tipo;
					if (tipo[0]!=' ' && tipo[1]=='\n'){bond_type << "     " << tipo[0];}
					if (tipo[0]!=' ' && tipo[1]!='\n')
					{
						if (tipo[1]==' '){bond_type << "     " << tipo [0];}
					}
					if (tipo[1]!=' ' && tipo[1]!='\n')
					{
						if (tipo[0]==' '){bond_type << "     " << tipo [1];}
						if (tipo[0]!=' '){bond_type << "    " << tipo[0] << tipo [1];}
					}
					if (tipo[2]!=' ' && tipo[2]!='\n')
					{
						if (tipo[1]==' '){bond_type << "     " << tipo [2];}
						if (tipo[1]!=' '){bond_type << "    " << tipo[1] << tipo [2];}
					}
					if (tipo[3]!=' ' && tipo[3]!='\n')
					{
						if (tipo[2]==' '){bond_type << "     " << tipo [3];}
						if (tipo[2]!=' '){bond_type << "    " << tipo[2] << tipo [3];}
					}
					if (tipo[4]!=' ' && tipo[4]!='\n')
					{
						if (tipo[3]==' '){bond_type << "     " << tipo [4];}
						if (tipo[3]!=' '){bond_type << "    " << tipo[3] << tipo [4];}
					}
					if (tipo[5]!=' ' && tipo[5]!='\n')
					{
						if (tipo[4]==' '){bond_type << "     " << tipo [5];}
						if (tipo[4]!=' '){bond_type << "    " << tipo[4] << tipo [5];}
					}
					if (tipo[6]!=' ' && tipo[6]!='\n')
					{
						if (tipo[5]==' '){bond_type << "     "  << tipo [6];}
						if (tipo[5]!=' '){bond_type << "    " << tipo[5] << tipo [6];}
					}
					bond_type << "\n";
					if (num_atomo_2<=9 && num_atomo_2!=0)
					{
						bonds << "     " << num_atomo_2 << "\n"; bond_type << "     " << num_atomo_2 << endl;
					}
					if (num_atomo_2<=99 && num_atomo_2>=10)
					{
						bonds << "    " << num_atomo_2 << "\n"; bond_type << "    " << num_atomo_2 << endl;
					}
					if (num_atomo_2<=999 && num_atomo_2>=100)
					{
						bonds << "   " << num_atomo_2 << "\n"; bond_type << "   " << num_atomo_2 << endl;
					}
					if (num_atomo_2<=9999 && num_atomo_2>=1000)
					{
						bonds << "  " << num_atomo_2 << "\n"; bond_type << "  " << num_atomo_2 << endl;
					}
					if (num_atomo_2<=99999 && num_atomo_2>=10000)
					{
						bonds << " " << num_atomo_2 << "\n";  bond_type << " " << num_atomo_2 << endl;
					}
					if (num_atomo_2<=999999 && num_atomo_2>=100000)
					{
						bonds << num_atomo_2 << "\n";  bond_type << num_atomo_2 << endl;
					}
					if (tipo[0]!=' ' && tipo[1]=='\n')
					{
						bond_type << "     " << tipo[0];
					}
					if (tipo[0]!=' ' && tipo[1]!='\n')
					{
					if (tipo[1]==' '){bond_type << "     " << tipo [0];}
					}
					if (tipo[1]!=' ' && tipo[1]!='\n'){
						if (tipo[0]==' ')
						{
							bond_type << "     "  << tipo [1];
						}
						if (tipo[0]!=' '){bond_type << "    " << tipo[0] << tipo [1];}
					}
					if (tipo[2]!=' ' && tipo[2]!='\n')
					{
						if (tipo[1]==' '){bond_type << "     " << tipo [2];}
						if (tipo[1]!=' '){bond_type << "    " << tipo[1] << tipo [2];}
					}
					if (tipo[3]!=' ' && tipo[3]!='\n')
					{
						if (tipo[2]==' '){bond_type << "     "  << tipo [3];}
						if (tipo[2]!=' '){bond_type << "    " << tipo[2] << tipo [3];}
					}
					if (tipo[4]!=' ' && tipo[4]!='\n')
					{
						if (tipo[3]==' '){bond_type << "     " << tipo [4];}
						if (tipo[3]!=' '){bond_type << "    " << tipo[3] << tipo [4];}
					}
					if (tipo[5]!=' ' && tipo[5]!='\n')
					{
						if (tipo[4]==' '){bond_type << "     " << tipo [5];}
						if (tipo[4]!=' '){bond_type << "    " << tipo[4] << tipo [5];}
					}
					bond_type << "\n";
				}//if (ordem[0]!='@' && atom_1[c1]!='@' && atom_2[c1]!='@')
				LAD.get(ordem[0]);
				if (!LAD.eof() && ordem[0]!='@' && ordem[0]!='\n'){goto point;}
				if (atomo_2[c1]==' '  && !LAD.eof())
				{
					line[53]=' ';	
					while (line[53]!='\n' && !LAD.eof())
					{
						LAD.get(line[53]);
					}
					LAD.get(line[53]);
				//	 if (line[53]=='\n'){while(!LAD.eof()){LAD.get(line[53]);}}
				}
				if (ordem[0]=='@' || ordem[c1]=='@' || atomo_1[c1]=='@' || atomo_2[c1]=='@' || ordem[c2]=='@' || atomo_1[c2]=='@' || atomo_2[c2]=='@' || line[53]=='@')
				{
					while (!LAD.eof()){LAD.get(line[0]);}
				}
				// cout << "ate aqui funcionou" << endl; //saida para ver se o programa conseguiu ler os caracteres ate este ponto
				//cout << atomo_2 << endl;   
				//else {cout << "erro" << ordem << endl << atomo_1 << endl << atomo_2 << endl << endl;}
			}
			//	 if (ordem[0]!=' '){while(!LAD.eof()){LAD.get(line[80]);}}
		} //while (pos!=18 && line[0]!='@' && !LAD.eof())
	}//while (!LAD.eof())
	LAD.close();
	bonds.close();	
	bond_type.close();
	//-------------------SEGUNDA-PARTE-----SAIDA------------------------
	ofstream root;
	root.open("L_ROOT_0.tmp");
	ifstream AD;
	AD.open("L_DIST.tmp");
	//--------variaveis do AD
	char AD_num_atomo[9];
	int num_AD_num_atomo;
	char AD_resto[46];
	//--------variaveis do AD
	//--------variaveis do bond
	char bond_1[8];
	char bond_2[7];
	int num_bond_1;
	int num_bond_2;
	//--------variaveis do bond
	//--------variaveis do L_root
	char espacinho[2];
	char L_root_1[8];
	char espacinho2[2];
	char L_root_atomo[4];
	char espacinho3[2];
	char L_root_resto[41];
	int num_L_root_1;
	//--------variaveis do L_root
	while (!AD.eof())
	{
		AD.get(AD_num_atomo,sizeof(AD_num_atomo),'\0');
		//cout << AD_num_atomo << endl;
		AD.get(AD_resto,sizeof(AD_resto),'\0');
		num_AD_num_atomo = atoi(AD_num_atomo);
		//cout << num_AD_num_atomo << endl << resto;
		ifstream bond;
		bond.open("bonds.tmp");
		while (!bond.eof())
		{
			bond.get(bond_1, sizeof(bond_1), '\0');
			//cout << bond_1;
			bond.get(bond_2, sizeof(bond_2), '\0');     
			//cout << bond_2;  
			num_bond_1 = atoi(bond_1);
			num_bond_2 = atoi(bond_2);
			if (num_bond_1==num_AD_num_atomo)
			{  
				ifstream L_root;
				L_root.open("ligand.tmp");
				while (!L_root.eof())
				{
					L_root.get(espacinho, sizeof(espacinho), '\0');
					L_root.get(L_root_1, sizeof(L_root_1), '\0');
					//cout << L_root_1 << endl;
					L_root.get(espacinho, sizeof(espacinho), '\0');
					L_root.get(L_root_atomo, sizeof(L_root_atomo), '\0');
					L_root.get(espacinho3, sizeof(espacinho3), '\0');
					L_root.get(L_root_resto, sizeof(L_root_resto), '\0');
					//cout << L_root_resto;
					num_L_root_1 = atoi(L_root_1);
					//cout << num_bond_1 << endl;
					//cout << num_bond_2 << endl;
					if (num_bond_2==num_L_root_1 && L_root_resto[34]!='H')
					{
						root << "R " << L_root_1 << L_root_atomo << "  " << L_root_resto; 
						//cout << num_bond_2 << endl << L_root_1 << L_root_resto << endl;
					}
				}// while (!L_root.eof())
				L_root.close();
			}//if (num_bond_1==num_AD_num_atomo)		
			if (num_bond_2==num_AD_num_atomo)
			{  
				ifstream L_root;
				L_root.open("ligand.tmp");
				while (!L_root.eof())
				{
					L_root.get(espacinho, sizeof(espacinho), '\0');
					L_root.get(L_root_1, sizeof(L_root_1), '\0');
					//cout << L_root_1 << endl;
					L_root.get(espacinho, sizeof(espacinho), '\0');
					L_root.get(L_root_atomo, sizeof(L_root_atomo), '\0');
					L_root.get(espacinho3, sizeof(espacinho3), '\0');
					L_root.get(L_root_resto, sizeof(L_root_resto), '\0');
					//cout << L_root_resto;
					num_L_root_1 = atoi(L_root_1);
					//cout << num_bond_1 << endl;
					//cout << num_bond_2 << endl;
					if (num_bond_1==num_L_root_1 && L_root_resto[34]!='H')
					{root << "R " << L_root_1 << " " << L_root_atomo << " " << L_root_resto;
					//cout << num_bond_2 << endl << L_root_1 << L_root_resto << endl;
					}
				}// while (!L_root.eof())
				L_root.close();
			}//if (num_bond_2==num_AD_num_atomo)
		}//while (!bonds.eof())
		bonds.close();
		root << "B" << AD_num_atomo << AD_resto;
	}//(!AD.eof())
	root.close();
	AD.close();
}//void raiz_ligante()

void P_centro_geom()
{
	//r1 = raiz1
	int pos = 0;
	int conta_r = 0;
	char tipo[2];
	char p_atomo[8];
	char p_nome_atomo[5];
	char p_aminoacido[4];
	char p_num_aminoacido[8];
	char px[9];
	char py[9];
	char pz[9];
	char p_enter[3];
	// raiz 1 
	char r1_atomo[9];
	char r1_nome_atomo[5];
	char r1_aminoacido[4];
	char r1_num_aminoacido[8];
	char r1x[9];
	char r1y[9];
	char r1z[9];
	char r1_enter[3];
	// raiz 2 
	char r2_atomo[9];
	char r2_nome_atomo[5];
	char r2_aminoacido[4];
	char r2_num_aminoacido[8];
	char r2x[9];
	char r2y[9];
	char r2z[9];
	char r2_enter[3];
	//--------converte caracteres para 2 raizes
	double num_r1x;
	double num_r1y;
	double num_r1z;
	double num_r2x;
	double num_r2y;
	double num_r2z;
	//--------converte caracteres para 3 raizes
	double num_r3x;
	double num_r3y;
	double num_r3z; 
	// raiz 3 
	char r3_atomo[9];
	char r3_nome_atomo[5];
	char r3_aminoacido[4];
	char r3_num_aminoacido[8];
	char r3x[9];
	char r3y[9];
	char r3z[9];
	char r3_enter[3];
	/* 
	// raiz 4 
	char r4_atomo[9];
	char r4_nome_atomo[5];
	char r4_aminoacido[4];
	char r4_num_aminoacido[8];
	char r4x[9];
	char r4y[9];
	char r4z[9];
	char r4_enter[3];
	*/
	//--------valores para calcular o centro geometrico (geometric center)
	double gcx2 = 0;
	double gcy2 = 0;
	double gcz2 = 0;
	double gcx2t = 0;
	double gcy2t = 0;
	double gcz2t = 0;
	double gcx3 = 0;
	double gcy3 = 0;
	double gcz3 = 0;
	double gcx3t = 0;
	double gcy3t = 0;
	double gcz3t = 0;
	ofstream p_root_f; //final
	
	
	p_root_f.open("P_ROOT.tmp"); 
	ifstream p_root;
	p_root.open("P_ROOT_0.tmp");
	if (!p_root)
	{
		cout << "Could not create temp file!" << endl;
	}

	while (!p_root.eof())
	{
		//pega uma linha 
		get:  p_root.get(tipo,sizeof(tipo),'\0');
		p_root.get(p_atomo,sizeof(p_atomo),'\0');
		p_root.get(p_nome_atomo,sizeof(p_nome_atomo),'\0');
		p_root.get(p_aminoacido,sizeof(p_aminoacido),'\0');
		p_root.get(p_num_aminoacido,sizeof(p_num_aminoacido),'\0');
		p_root.get(px,sizeof(px),'\0');
		p_root.get(py,sizeof(py),'\0');
		p_root.get(pz,sizeof(pz),'\0');
		p_root.get(p_enter,sizeof(p_enter),'\0');
		// se nao for raiz e tiver so uma raiz
		if (tipo[0]=='B' && conta_r==1)
		{
			conta_r = 0;
			p_root_f << "R" << r1_atomo << r1_nome_atomo << r1_aminoacido << r1_num_aminoacido << r1x << r1y << r1z << r1_enter;
			p_root_f << "B" << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
			goto get;
		}
		//se a primeira linha for uma raiz  
		if (tipo[0]=='R' && conta_r==0)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r1_atomo[pos] = p_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r1_nome_atomo[pos] = p_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=5)
			{r1_aminoacido[pos] = p_aminoacido[pos]; pos++;}
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r1_num_aminoacido[pos] = p_num_aminoacido[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=10)
			{r1x[pos] = px[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=10)
			{r1y[pos] = py[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=10)
			{r1z[pos] = pz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=5)
			{r1_enter[pos] = p_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r1_atomo << "*" << r1_nome_atomo << "*" << r1_aminoacido << "*" << r1_num_aminoacido << "*" << r1x << "*" << r1y << "*" << r1z << "*" << r1_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		//se a primeira linha for uma raiz  
		if (tipo[0]=='R' && conta_r==1)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r2_atomo[pos] = p_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r2_nome_atomo[pos] = p_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=5)
			{r2_aminoacido[pos] = p_aminoacido[pos]; pos++;}
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r2_num_aminoacido[pos] = p_num_aminoacido[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=10)
			{r2x[pos] = px[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=10)
			{r2y[pos] = py[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=10)
			{r2z[pos] = pz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=5)
			{r2_enter[pos] = p_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r2_atomo << "*" << r2_nome_atomo << "*" << r2_aminoacido << "*" << r2_num_aminoacido << "*" << r2x << "*" << r2y << "*" << r2z << "*" << r2_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		if (tipo[0]=='R' && conta_r==2)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r3_atomo[pos] = p_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r3_nome_atomo[pos] = p_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=5)
			{r3_aminoacido[pos] = p_aminoacido[pos]; pos++;}
			//------copia do numero do atomo
			pos = 0;
			while (pos!=9)
			{r3_num_aminoacido[pos] = p_num_aminoacido[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=10)
			{r3x[pos] = px[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=10)
			{r3y[pos] = py[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=10)
			{r3z[pos] = pz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=5)
			{r3_enter[pos] = p_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r3_atomo << "*" << r3_nome_atomo << "*" << r3_aminoacido << "*" << r3_num_aminoacido << "*" << r3x << "*" << r3y << "*" << r3z << "*" << r3_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		// se tiver 2 ou 3 raizes
		if (tipo[0]=='B' && (conta_r==2 || conta_r==3))
		{
				//elimina erros dos vetores da primeira raiz
				if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; }
				if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; }
				if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; }
				if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; }
				if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; }
				if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; }
				if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; }
				if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; }
				if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; }
				if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; }
				if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; }
				if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; }
				//elimina erros dos vetores da segunda raiz
				if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; }
				if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; }
				if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; }
				if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; }
				if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; }
				if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; }
				if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; }
				if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; }
				if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; }
				if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; }
				if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; }
				if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; }
				num_r1x = atof(r1x);
				num_r1y = atof(r1y);
				num_r1z = atof(r1z);
				//cout << num_r1x << endl << num_r1y << endl << num_r1z << endl;
				num_r2x = atof(r2x);
				num_r2y = atof(r2y);
				num_r2z = atof(r2z);   		 
				//cout << num_r2x << endl << num_r2y << endl << num_r2z << endl;
				if (conta_r==2)
				{   	
					gcx2t = num_r1x + num_r2x;
					gcx2 = gcx2t/2;
					gcy2t = num_r1y + num_r2y;
					gcy2 = gcy2t/2;
					gcz2t = num_r1z + num_r2z;
					gcz2 = gcz2t/2;	
					//cout << "resultado(2): " << gcx3 << " / " << gcy3 << " / " << gcz3 << endl;
					p_root_f << "R   " << "geom. center(2)   ";
					char letrasx[11];
					char letrasy[11];
					char letrasz[11];
					letrasx[0] = ' ';
					letrasx[1] = ' ';
					letrasx[2] = ' ';
					letrasx[3] = ' ';
					letrasx[4] = ' ';
					letrasx[5] = ' ';
					letrasx[6] = ' ';
					letrasx[7] = ' ';
					letrasx[8] = ' ';
					letrasx[9] = ' ';
					letrasy[0] = ' ';
					letrasy[1] = ' ';
					letrasy[2] = ' ';
					letrasy[3] = ' ';
					letrasy[4] = ' ';
					letrasy[5] = ' ';
					letrasy[6] = ' ';
					letrasy[7] = ' ';
					letrasy[8] = ' ';
					letrasy[9] = ' ';
					letrasz[0] = ' ';
					letrasz[1] = ' ';
					letrasz[2] = ' ';
					letrasz[3] = ' ';
					letrasz[4] = ' ';
					letrasz[5] = ' ';
					letrasz[6] = ' ';
					letrasz[7] = ' ';
					letrasz[8] = ' ';
					letrasz[9] = ' ';
					if (gcx2>=-99.9999 && gcx2<=999.9999)
					{sprintf(letrasx, "%.4f", gcx2);}
					if (gcx2>=-999.9999 && gcx2<=-100.0000)
					{sprintf(letrasx, "%.3f", gcx2);}
					if (gcx2>=-99.9999 && gcy2<=999.9999)
					{sprintf(letrasy, "%.4f", gcy2);}
					if (gcy2>=-999.9999 && gcy2<=-100.0000)
					{sprintf(letrasy, "%.3f", gcy2);}
					if (gcz2>=-99.9999 && gcz2<=999.9999)
					{sprintf(letrasz, "%.4f", gcz2);}
					if (gcz2>=-999.9999 && gcz2<=-100.0000)
					{sprintf(letrasz, "%.3f", gcz2);}
					//---------------X
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					//---------------Y
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					//---------------Z
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 				    		 
					//cout << letrasx << endl;
					//cout << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
					// p_root_f << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7]; 
					p_root_f << letrasx << letrasy << letrasz;  		 		 
					p_root_f << " " << endl;
					conta_r = 0;
				}
				if (conta_r==3)
				{
					//elimina erros dos vetores da segunda raiz
					if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; }
					if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; }
					if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; }
					if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; }
					if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; }
					if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; }
					if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; }
					if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; }
					if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; }
					if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; }
					if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; }
					if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; }
					
					################### parei aqui ####################################################################
					num_r3x = atof(r3x);
					num_r3y = atof(r3y);
					num_r3z = atof(r3z);   		 
					gcx3t = (num_r1x+num_r2x+num_r3x);
					gcx3 = (gcx3t/3);
					gcy3t = (num_r1y+num_r2y+num_r3y);
					gcy3 = (gcy3t/3);
					gcz3t = (num_r1z+num_r2z+num_r3z);
					gcz3 = (gcz3t/3);
					/*
					cout << num_r3x << endl << num_r3y << endl << num_r3z << endl;
					cout << endl << "soma temp" << endl;
					cout << gcx3t << endl << gcy3t << endl << gcz3t << endl << endl;
					*/
					//gcx3 = -50.2211; gcy3 = -100.2234; gcz3 = -0.3345;	 
					p_root_f << "R   " << "geom. center(3)   ";
					char letrasx[11];
					char letrasy[11];
					char letrasz[11];
					letrasx[0] = ' ';
					letrasx[1] = ' ';
					letrasx[2] = ' ';
					letrasx[3] = ' ';
					letrasx[4] = ' ';
					letrasx[5] = ' ';
					letrasx[6] = ' ';
					letrasx[7] = ' ';
					letrasx[8] = ' ';
					letrasx[9] = ' ';
					letrasy[0] = ' ';
					letrasy[1] = ' ';
					letrasy[2] = ' ';
					letrasy[3] = ' ';
					letrasy[4] = ' ';
					letrasy[5] = ' ';
					letrasy[6] = ' ';
					letrasy[7] = ' ';
					letrasy[8] = ' ';
					letrasy[9] = ' ';
					letrasz[0] = ' ';
					letrasz[1] = ' ';
					letrasz[2] = ' ';
					letrasz[3] = ' ';
					letrasz[4] = ' ';
					letrasz[5] = ' ';
					letrasz[6] = ' ';
					letrasz[7] = ' ';
					letrasz[8] = ' ';
					letrasz[9] = ' ';
					if (gcx3>=-99.9999 && gcx3<=999.9999)
					{sprintf(letrasx, "%.4f", gcx3);}
					if (gcx3>=-999.9999 && gcx3<=-100.0000)
					{sprintf(letrasx, "%.3f", gcx3);}
					if (gcx3>=-99.9999 && gcy3<=999.9999)
					{sprintf(letrasy, "%.4f", gcy3);}
					if (gcy3>=-999.9999 && gcy3<=-100.0000)
					{sprintf(letrasy, "%.3f", gcy3);}
					if (gcz3>=-99.9999 && gcz3<=999.9999)
					{sprintf(letrasz, "%.4f", gcz3);}
					if (gcz3>=-999.9999 && gcz3<=-100.0000)
					{sprintf(letrasz, "%.3f", gcz3);}
					//---------------X
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					if (letrasx[8]==' ')
					{
						letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
					} 
					//---------------Y
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					if (letrasy[8]==' ')
					{
						letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
					} 
					//---------------Z
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 
					if (letrasz[8]==' ')
					{
						letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
					} 				    		 
					//cout << letrasx << endl;
					//cout << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
					// p_root_f << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7]; 
					p_root_f << letrasx << letrasy << letrasz;  		 		 
					p_root_f << " " << endl;
					//  p_root_f << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
					conta_r = 0;
				}//if (conta_r==3)	 
				p_root_f << tipo << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
		}//if (conta_r==2 || conta_r==3)
		if (tipo[0]!='R' && tipo[0]!='B' && tipo[0]!='\0')
		{
		cout << "ERRO! |" << tipo <<"|" << endl;  	 			   
		}
		//cout << tipo << p_atomo << p_nome_atomo << p_aminoacido << p_num_aminoacido << px << py << pz << p_enter;
	}//while (!p_root.eof())
	p_root.close(); 
	p_root_f.close();
}

//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
void L_centro_geom()
{
	//r1 = raiz1
	int pos = 0;
	int conta_r = 0;
	char tipo[2];
	char l_atomo[9];
	char l_nome_atomo[5];
	char l_resto[6];
	char lx[11];
	char ly[11];
	char lz[11];
	char l_enter[7];
	// raiz 1 
	char r1_atomo[9];
	char r1_nome_atomo[5];
	char r1_resto[6];
	char r1x[11];
	char r1y[11];
	char r1z[11];
	char r1_enter[7];
	// raiz 2 
	char r2_atomo[9];
	char r2_nome_atomo[5];
	char r2_resto[6];
	char r2x[11];
	char r2y[11];
	char r2z[11];
	char r2_enter[7];
	//--------converte caracteres para 2 raizes
	double num_r1x;
	double num_r1y;
	double num_r1z;
	double num_r2x;
	double num_r2y;
	double num_r2z;
	//--------converte caracteres para 3 raizes
	double num_r3x;
	double num_r3y;
	double num_r3z; 
	// raiz 3 
	char r3_atomo[9];
	char r3_nome_atomo[5];
	char r3_resto[6]; 
	char r3x[11];
	char r3y[11];
	char r3z[11];
	char r3_enter[7];
	/* 
	// raiz 4 
	char r4_atomo[9];
	char r4_nome_atomo[5];
	char r4_aminoacido[4];
	char r4_num_aminoacido[8];
	char r4x[9];
	char r4y[9];
	char r4z[9];
	char r4_enter[3];
	*/
	//--------valores para calcular o centro geometrico (geometric center)
	double gcx2 = 0;
	double gcy2 = 0;
	double gcz2 = 0;
	double gcx2t = 0;
	double gcy2t = 0;
	double gcz2t = 0;
	double gcx3 = 0;
	double gcy3 = 0;
	double gcz3 = 0;
	double gcx3t = 0;
	double gcy3t = 0;
	double gcz3t = 0;
	ofstream l_root_f; //final
	l_root_f.open("L_ROOT.tmp"); 
	ifstream l_root;
	l_root.open("L_ROOT_0.tmp");
	if (!l_root)
	{cout << "Could not create temp file!" << endl;}
	while (!l_root.eof())
	{
		//pega uma linha 
		get:  l_root.get(tipo,sizeof(tipo),'\0');
		l_root.get(l_atomo,sizeof(l_atomo),'\0');
		l_root.get(l_nome_atomo,sizeof(l_nome_atomo),'\0');
		l_root.get(l_resto,sizeof(l_resto),'\0');
		l_root.get(lx,sizeof(lx),'\0');
		l_root.get(ly,sizeof(ly),'\0');
		l_root.get(lz,sizeof(lz),'\0');
		l_root.get(l_enter,sizeof(l_enter),'\0');
		//cout << "*" << tipo << "*" << l_atomo << "*" << l_nome_atomo << "*" << l_resto << "*" << lx << "*" << ly << "*" << lz << "*" << l_enter;
		// se nao for raiz e tiver so uma raiz
		if (tipo[0]=='B' && conta_r==1)
		{
			conta_r = 0;
			l_root_f << "R" << r1_atomo << r1_nome_atomo << r1_resto << r1x << r1y << r1z << " " << endl;
			l_root_f << "B" << l_atomo << l_nome_atomo << l_resto  << lx << ly << lz << " " << endl;

			goto get;
		}
		//se a primeira linha for uma raiz  
		if (tipo[0]=='R' && conta_r==0)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=10)
			{r1_atomo[pos] = l_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r1_nome_atomo[pos] = l_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=7)
			{r1_resto[pos] = l_resto[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=12)
			{r1x[pos] = lx[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=12)
			{r1y[pos] = ly[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=12)
			{r1z[pos] = lz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=8)
			{r1_enter[pos] = l_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r1_atomo << "*" << r1_nome_atomo << "*" << r1_aminoacido << "*" << r1_num_aminoacido << "*" << r1x << "*" << r1y << "*" << r1z << "*" << r1_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		//se a primeira linha for uma raiz  
		if (tipo[0]=='R' && conta_r==1)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=10)
			{r2_atomo[pos] = l_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r2_nome_atomo[pos] = l_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=7)
			{r2_resto[pos] = l_resto[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=12)
			{r2x[pos] = lx[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=12)
			{r2y[pos] = ly[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=12)
			{r2z[pos] = lz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=8)
			{r2_enter[pos] = l_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r2_atomo << "*" << r2_nome_atomo << "*" << r2_aminoacido << "*" << r2_num_aminoacido << "*" << r2x << "*" << r2y << "*" << r2z << "*" << r2_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		if (tipo[0]=='R' && conta_r==2)
		{
			conta_r++;
			//------copia do numero do atomo
			pos = 0;
			while (pos!=10)
			{r3_atomo[pos] = l_atomo[pos]; pos++;}
			//------copia do nome do atomo
			pos = 0;
			while (pos!=6)
			{r3_nome_atomo[pos] = l_nome_atomo[pos]; pos++;}
			//------copia do nome do aa
			pos = 0;
			while (pos!=7)
			{r3_resto[pos] = l_resto[pos]; pos++;}
			//------copia posicao x
			pos = 0;
			while (pos!=12)
			{r3x[pos] = lx[pos]; pos++;}
			//------copia posicao y
			pos = 0;
			while (pos!=12)
			{r3y[pos] = ly[pos]; pos++;}
			//------copia posicao z
			pos = 0;
			while (pos!=12)
			{r3z[pos] = lz[pos]; pos++;}
			//------copia final
			pos = 0;
			while (pos!=8)
			{r3_enter[pos] = l_enter[pos]; pos++;}
			pos = 0;
			//cout << "*" << r3_atomo << "*" << r3_nome_atomo << "*" << r3_aminoacido << "*" << r3_num_aminoacido << "*" << r3x << "*" << r3y << "*" << r3z << "*" << r3_enter;
			//cout << endl << conta_r << endl;  
			goto get;
		}//if (tipo[0]=='R' && conta_r==0)
		// se tiver 2 ou 3 raizes
		if (tipo[0]=='B' && (conta_r==2 || conta_r==3))
		{
			//elimina erros dos vetores da primeira raiz
			if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; r1x[10]==r1x[11]; }
			if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; r1x[10]==r1x[11]; }
			if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; r1x[10]==r1x[11]; }
			if (r1x[0]==' '){r1x[0]==r1x[1]; r1x[1]==r1x[2]; r1x[2]==r1x[3]; r1x[3]==r1x[4]; r1x[4]==r1x[5]; r1x[5]==r1x[6]; r1x[6]==r1x[7]; r1x[8]==r1x[9]; r1x[9]==r1x[10]; r1x[10]==r1x[11]; }
			if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; r1y[10]==r1y[11]; }
			if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; r1y[10]==r1y[11]; }
			if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; r1y[10]==r1y[11]; }
			if (r1y[0]==' '){r1y[0]==r1y[1]; r1y[1]==r1y[2]; r1y[2]==r1y[3]; r1y[3]==r1y[4]; r1y[4]==r1y[5]; r1y[5]==r1y[6]; r1y[6]==r1y[7]; r1y[8]==r1y[9]; r1y[9]==r1y[10]; r1y[10]==r1y[11]; }
			if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; r1z[10]==r1z[11]; }
			if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; r1z[10]==r1z[11]; }
			if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; r1z[10]==r1z[11]; }
			if (r1z[0]==' '){r1z[0]==r1z[1]; r1z[1]==r1z[2]; r1z[2]==r1z[3]; r1z[3]==r1z[4]; r1z[4]==r1z[5]; r1z[5]==r1z[6]; r1z[6]==r1z[7]; r1z[8]==r1z[9]; r1z[9]==r1z[10]; r1z[10]==r1z[11]; }
			//elimina erros dos vetores da segunda raiz
			if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; r2x[10]==r2x[11]; }
			if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; r2x[10]==r2x[11]; }
			if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; r2x[10]==r2x[11]; }
			if (r2x[0]==' '){r2x[0]==r2x[1]; r2x[1]==r2x[2]; r2x[2]==r2x[3]; r2x[3]==r2x[4]; r2x[4]==r2x[5]; r2x[5]==r2x[6]; r2x[6]==r2x[7]; r2x[8]==r2x[9]; r2x[9]==r2x[10]; r2x[10]==r2x[11]; }
			if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; r2y[10]==r2y[11]; }
			if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; r2y[10]==r2y[11]; }
			if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; r2y[10]==r2y[11]; }
			if (r2y[0]==' '){r2y[0]==r2y[1]; r2y[1]==r2y[2]; r2y[2]==r2y[3]; r2y[3]==r2y[4]; r2y[4]==r2y[5]; r2y[5]==r2y[6]; r2y[6]==r2y[7]; r2y[8]==r2y[9]; r2y[9]==r2y[10]; r2y[10]==r2y[11]; }
			if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; r2z[10]==r2z[11]; }
			if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; r2z[10]==r2z[11]; }
			if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; r2z[10]==r2z[11]; }
			if (r2z[0]==' '){r2z[0]==r2z[1]; r2z[1]==r2z[2]; r2z[2]==r2z[3]; r2z[3]==r2z[4]; r2z[4]==r2z[5]; r2z[5]==r2z[6]; r2z[6]==r2z[7]; r2z[8]==r2z[9]; r2z[9]==r2z[10]; r2z[10]==r2z[11]; }
			num_r1x = atof(r1x);
			num_r1y = atof(r1y);
			num_r1z = atof(r1z);
			//cout << num_r1x << endl << num_r1y << endl << num_r1z << endl;
			num_r2x = atof(r2x);
			num_r2y = atof(r2y);
			num_r2z = atof(r2z);   		 
			//cout << num_r2x << endl << num_r2y << endl << num_r2z << endl;
			if (conta_r==2)
			{   	
				gcx2t = num_r1x + num_r2x;
				gcx2 = gcx2t/2;
				gcy2t = num_r1y + num_r2y;
				gcy2 = gcy2t/2;
				gcz2t = num_r1z + num_r2z;
				gcz2 = gcz2t/2;	
				//cout << "resultado(2): " << gcx3 << " / " << gcy3 << " / " << gcz3 << endl;
				l_root_f << "R " << "geom. center(2)  ";
				char letrasx[11];
				char letrasy[11];
				char letrasz[11];
				letrasx[0] = ' ';
				letrasx[1] = ' ';
				letrasx[2] = ' ';
				letrasx[3] = ' ';
				letrasx[4] = ' ';
				letrasx[5] = ' ';
				letrasx[6] = ' ';
				letrasx[7] = ' ';
				letrasx[8] = ' ';
				letrasx[9] = ' ';
				letrasy[0] = ' ';
				letrasy[1] = ' ';
				letrasy[2] = ' ';
				letrasy[3] = ' ';
				letrasy[4] = ' ';
				letrasy[5] = ' ';
				letrasy[6] = ' ';
				letrasy[7] = ' ';
				letrasy[8] = ' ';
				letrasy[9] = ' ';
				letrasz[0] = ' ';
				letrasz[1] = ' ';
				letrasz[2] = ' ';
				letrasz[3] = ' ';
				letrasz[4] = ' ';
				letrasz[5] = ' ';
				letrasz[6] = ' ';
				letrasz[7] = ' ';
				letrasz[8] = ' ';
				letrasz[9] = ' ';
				if (gcx2>=-99.9999 && gcx2<=999.9999)
				{sprintf(letrasx, "%.4f", gcx2);}
				if (gcx2>=-999.9999 && gcx2<=-100.0000)
				{sprintf(letrasx, "%.3f", gcx2);}
				if (gcx2>=-99.9999 && gcy2<=999.9999)
				{sprintf(letrasy, "%.4f", gcy2);}
				if (gcy2>=-999.9999 && gcy2<=-100.0000)
				{sprintf(letrasy, "%.3f", gcy2);}
				if (gcz2>=-99.9999 && gcz2<=999.9999)
				{sprintf(letrasz, "%.4f", gcz2);}
				if (gcz2>=-999.9999 && gcz2<=-100.0000)
				{sprintf(letrasz, "%.3f", gcz2);}
				//---------------X
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				//---------------Y
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
				} 
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
				} 
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
				} 
				//---------------Z
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 				    		 
				//cout << letrasx << endl;
				//cout << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
				// p_root_f << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7]; 
				l_root_f << letrasx << "  " << letrasy << "  " << letrasz;  		 		 
				l_root_f << "  " << endl;
				conta_r = 0;
			}
			if (conta_r==3)
			{
				//elimina erros dos vetores da segunda raiz
				if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; r3x[10]==r3x[11]; }
				if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; r3x[10]==r3x[11]; }
				if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; r3x[10]==r3x[11]; }
				if (r3x[0]==' '){r3x[0]==r3x[1]; r3x[1]==r3x[2]; r3x[2]==r3x[3]; r3x[3]==r3x[4]; r3x[4]==r3x[5]; r3x[5]==r3x[6]; r3x[6]==r3x[7]; r3x[8]==r3x[9]; r3x[9]==r3x[10]; r3x[10]==r3x[11]; }
				if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; r3y[10]==r3y[11]; }
				if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; r3y[10]==r3y[11]; }
				if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; r3y[10]==r3y[11]; }
				if (r3y[0]==' '){r3y[0]==r3y[1]; r3y[1]==r3y[2]; r3y[2]==r3y[3]; r3y[3]==r3y[4]; r3y[4]==r3y[5]; r3y[5]==r3y[6]; r3y[6]==r3y[7]; r3y[8]==r3y[9]; r3y[9]==r3y[10]; r3y[10]==r3y[11]; }
				if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; r3z[10]==r3z[11]; }
				if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; r3z[10]==r3z[11]; }
				if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; r3z[10]==r3z[11]; }
				if (r3z[0]==' '){r3z[0]==r3z[1]; r3z[1]==r3z[2]; r3z[2]==r3z[3]; r3z[3]==r3z[4]; r3z[4]==r3z[5]; r3z[5]==r3z[6]; r3z[6]==r3z[7]; r3z[8]==r3z[9]; r3z[9]==r3z[10]; r3z[10]==r3z[11]; }
				num_r3x = atof(r3x);
				num_r3y = atof(r3y);
				num_r3z = atof(r3z);   		 
				gcx3t = (num_r1x+num_r2x+num_r3x);
				gcx3 = (gcx3t/3);
				gcy3t = (num_r1y+num_r2y+num_r3y);
				gcy3 = (gcy3t/3);
				gcz3t = (num_r1z+num_r2z+num_r3z);
				gcz3 = (gcz3t/3);
				//cout << num_r3x << endl << num_r3y << endl << num_r3z << endl;
				//cout << endl << "soma temp" << endl;
				//cout << gcx3t << endl << gcy3t << endl << gcz3t << endl << endl;
				//gcx3 = -50.2211; gcy3 = -100.2234; gcz3 = -0.3345;	 
				l_root_f << "R " << "geom. center(3)  ";
				char letrasx[11];
				char letrasy[11];
				char letrasz[11];
				letrasx[0] = ' ';
				letrasx[1] = ' ';
				letrasx[2] = ' ';
				letrasx[3] = ' ';
				letrasx[4] = ' ';
				letrasx[5] = ' ';
				letrasx[6] = ' ';
				letrasx[7] = ' ';
				letrasx[8] = ' ';
				letrasx[9] = ' ';
				letrasy[0] = ' ';
				letrasy[1] = ' ';
				letrasy[2] = ' ';
				letrasy[3] = ' ';
				letrasy[4] = ' ';
				letrasy[5] = ' ';
				letrasy[6] = ' ';
				letrasy[7] = ' ';
				letrasy[8] = ' ';
				letrasy[9] = ' ';
				letrasz[0] = ' ';
				letrasz[1] = ' ';
				letrasz[2] = ' ';
				letrasz[3] = ' ';
				letrasz[4] = ' ';
				letrasz[5] = ' ';
				letrasz[6] = ' ';
				letrasz[7] = ' ';
				letrasz[8] = ' ';
				letrasz[9] = ' ';
				if (gcx3>=-99.9999 && gcx3<=999.9999)
				{sprintf(letrasx, "%.4f", gcx3);}
				if (gcx3>=-999.9999 && gcx3<=-100.0000)
				{sprintf(letrasx, "%.3f", gcx3);}
				if (gcx3>=-99.9999 && gcy3<=999.9999)
				{sprintf(letrasy, "%.4f", gcy3);}
				if (gcy3>=-999.9999 && gcy3<=-100.0000)
				{sprintf(letrasy, "%.3f", gcy3);}
				if (gcz3>=-99.9999 && gcz3<=999.9999)
				{sprintf(letrasz, "%.4f", gcz3);}
				if (gcz3>=-999.9999 && gcz3<=-100.0000)
				{sprintf(letrasz, "%.3f", gcz3);}
				//---------------X
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				if (letrasx[8]==' ')
				{
					letrasx[8]=letrasx[7]; letrasx[7]=letrasx[6]; letrasx[6]=letrasx[5]; letrasx[5]=letrasx[4]; letrasx[4]=letrasx[3]; letrasx[3]=letrasx[2]; letrasx[2]=letrasx[1]; letrasx[1]=letrasx[0]; letrasx[0]=' ';			 
				} 
				//---------------Y
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
				} 
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';			 
				} 
				if (letrasy[8]==' ')
				{
					letrasy[8]=letrasy[7]; letrasy[7]=letrasy[6]; letrasy[6]=letrasy[5]; letrasy[5]=letrasy[4]; letrasy[4]=letrasy[3]; letrasy[3]=letrasy[2]; letrasy[2]=letrasy[1]; letrasy[1]=letrasy[0]; letrasy[0]=' ';
				} 
				//---------------Z
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 
				if (letrasz[8]==' ')
				{
					letrasz[8]=letrasz[7]; letrasz[7]=letrasz[6]; letrasz[6]=letrasz[5]; letrasz[5]=letrasz[4]; letrasz[4]=letrasz[3]; letrasz[3]=letrasz[2]; letrasz[2]=letrasz[1]; letrasz[1]=letrasz[0]; letrasz[0]=' ';			 
				} 				    		 
				//cout << letrasx << endl;
				//cout << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
				// p_root_f << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7]; 
				l_root_f << letrasx << "  " << letrasy << "  " << letrasz;  		 		 
				l_root_f << "  " << endl;
				//  p_root_f << l_atomo << l_nome_atomo << p_aminoacido << p_num_aminoacido << lx << ly << lz << l_enter;
				conta_r = 0;
			}//if (conta_r==3)	 
			l_root_f << tipo << l_atomo << l_nome_atomo << l_resto << lx << ly << lz << " " << endl;
		}//if (conta_r==2 || conta_r==3)
		if (tipo[0]!='R' && tipo[0]!='B' && tipo[0]!='\0')
		{
			cout << "ERROR! Can't find ligand/protein atoms root! |" << tipo <<"|" << endl;  	 			   
		}
		//cout << tipo << l_atomo << l_nome_atomo << p_aminoacido << p_num_aminoacido << lx << ly << lz << l_enter;
	}//while (!p_root.eof())
	l_root.close(); 
	l_root_f.close();
}

//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
void angulos()
{
	int bond_count = 0;
	//--------variaveis para carregar a proteina
	//raiz
	char rp_tipo[2];
	char rp_atomo[8];
	char rp_nome_atomo[5];
	char rp_aminoacido[4];
	char rp_num_aminoacido[8];
	char c1_rpx[9];
	char c1_rpy[9];
	char c1_rpz[9];
	char rpx[9];
	char rpy[9];
	char rpz[9];
	double num_rpx;
	double num_rpy;
	double num_rpz;
	char rp_enter[3];

	//aceitador/doador
	char bp_tipo[2];
	char bp_atomo[8];
	char bp_nome_atomo[5];
	char bp_aminoacido[4];
	char bp_num_aminoacido[8];
	char c1_bpx[9];
	char c1_bpy[9];
	char c1_bpz[9];
	char bpx[9];
	char bpy[9];
	char bpz[9];
	double num_bpx;
	double num_bpy;
	double num_bpz;
	char bp_enter[3];

	//--------------------------- 
	//---------variaveis para carregar o ligante
	char rl_tipo[2];
	char rl_atomo[9];
	char rl_nome_atomo[5];
	char rl_resto[6];
	char c1_rlx[11];
	char c1_rly[11];
	char c1_rlz[11];
	char rlx[11];
	char rly[11];
	char rlz[11];
	double num_rlx;
	double num_rly;
	double num_rlz;
	char rl_enter[3];
	char bl_tipo[2];
	char bl_atomo[9];
	char bl_nome_atomo[5];
	char bl_resto[6];
	char c1_blx[11];
	char c1_bly[11];
	char c1_blz[11];
	char blx[11];
	char bly[11];
	char blz[11];
	double num_blx;
	double num_bly;
	double num_blz;
	char bl_enter[3];

	//------saida-arquivos-temporarios com ligante e proteina (atomos mais proximos)
	ofstream p_result;
	p_result.open("p_result.tmp");
	ofstream l_result;
	l_result.open("l_result.tmp");
	ofstream distances;
	distances.open("dist_result.tmp");

	//-----------------------
	ifstream proteina;
	ifstream ligante;
	proteina.open("P_ROOT.tmp");
	ligante.open("L_ROOT.tmp");
	while (!proteina.eof() || !ligante.eof()){
		proteina.get(rp_tipo,sizeof(rp_tipo),'\0');
		proteina.get(rp_atomo,sizeof(rp_atomo),'\0');
		proteina.get(rp_nome_atomo,sizeof(rp_nome_atomo),'\0');
		proteina.get(rp_aminoacido,sizeof(rp_aminoacido),'\0');
		proteina.get(rp_num_aminoacido,sizeof(rp_num_aminoacido),'\0');
		proteina.get(rpx,sizeof(rpx),'\0');
		proteina.get(rpy,sizeof(rpy),'\0');
		proteina.get(rpz,sizeof(rpz),'\0');
		proteina.get(rp_enter,sizeof(rp_enter),'\0');
        
		proteina.get(bp_tipo,sizeof(bp_tipo),'\0');
		proteina.get(bp_atomo,sizeof(bp_atomo),'\0');
		proteina.get(bp_nome_atomo,sizeof(bp_nome_atomo),'\0');
		proteina.get(bp_aminoacido,sizeof(bp_aminoacido),'\0');
		proteina.get(bp_num_aminoacido,sizeof(bp_num_aminoacido),'\0');
		proteina.get(bpx,sizeof(bpx),'\0');
		proteina.get(bpy,sizeof(bpy),'\0');
		proteina.get(bpz,sizeof(bpz),'\0');
		proteina.get(bp_enter,sizeof(bp_enter),'\0');
		//	cout << "*" << p_tipo << "*" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << px << "*" << py << "*" << pz << "*" << p_enter;					
		ligante.get(rl_tipo,sizeof(rl_tipo),'\0');
		ligante.get(rl_atomo,sizeof(rl_atomo),'\0');
		ligante.get(rl_nome_atomo,sizeof(rl_nome_atomo),'\0');
		ligante.get(rl_resto,sizeof(rl_resto),'\0');
		ligante.get(rlx,sizeof(rlx),'\0');
		ligante.get(rly,sizeof(rly),'\0');
		ligante.get(rlz,sizeof(rlz),'\0');
		ligante.get(rl_enter,sizeof(rl_enter),'\0');
		ligante.get(bl_tipo,sizeof(bl_tipo),'\0');
		ligante.get(bl_atomo,sizeof(bl_atomo),'\0');
		ligante.get(bl_nome_atomo,sizeof(bl_nome_atomo),'\0');
		ligante.get(bl_resto,sizeof(bl_resto),'\0');
		ligante.get(blx,sizeof(blx),'\0');
		ligante.get(bly,sizeof(bly),'\0');
		ligante.get(blz,sizeof(blz),'\0');
		ligante.get(bl_enter,sizeof(bl_enter),'\0');
		// cout << "*" << l_tipo << "*" << l_atomo << "*" << l_nome_atomo << "*" << l_resto << "*" << lx << "*" << ly << "*" << lz << "*" << l_enter;					
		//--------antes de comecar a alterar o xyz copia string para outra
		// raiz proteina
		strcpy(c1_rpx, rpx);
		strcpy(c1_rpy, rpy);
		strcpy(c1_rpz, rpz);
		// raiz ligante
		strcpy(c1_rlx, rlx);
		strcpy(c1_rly, rly);
		strcpy(c1_rlz, rlz);
		// aceitador / doador proteina
		strcpy(c1_bpx, bpx);
		strcpy(c1_bpy, bpy);
		strcpy(c1_bpz, bpz);
		// aceitador / doador ligante
		strcpy(c1_blx, blx);
		strcpy(c1_bly, bly);
		strcpy(c1_blz, blz);
		//----------------------
		//elimina erros dos vetores da raiz da proteina
		if (rpx[0]==' '){rpx[0]==rpx[1]; rpx[1]==rpx[2]; rpx[2]==rpx[3]; rpx[3]==rpx[4]; rpx[4]==rpx[5]; rpx[5]==rpx[6]; rpx[6]==rpx[7]; rpx[8]==rpx[9]; rpx[9]==rpx[10]; rpx[10]==rpx[11]; }
		if (rpx[0]==' '){rpx[0]==rpx[1]; rpx[1]==rpx[2]; rpx[2]==rpx[3]; rpx[3]==rpx[4]; rpx[4]==rpx[5]; rpx[5]==rpx[6]; rpx[6]==rpx[7]; rpx[8]==rpx[9]; rpx[9]==rpx[10]; rpx[10]==rpx[11]; }
		if (rpx[0]==' '){rpx[0]==rpx[1]; rpx[1]==rpx[2]; rpx[2]==rpx[3]; rpx[3]==rpx[4]; rpx[4]==rpx[5]; rpx[5]==rpx[6]; rpx[6]==rpx[7]; rpx[8]==rpx[9]; rpx[9]==rpx[10]; rpx[10]==rpx[11]; }
		if (rpx[0]==' '){rpx[0]==rpx[1]; rpx[1]==rpx[2]; rpx[2]==rpx[3]; rpx[3]==rpx[4]; rpx[4]==rpx[5]; rpx[5]==rpx[6]; rpx[6]==rpx[7]; rpx[8]==rpx[9]; rpx[9]==rpx[10]; rpx[10]==rpx[11]; }
		if (rpy[0]==' '){rpy[0]==rpy[1]; rpy[1]==rpy[2]; rpy[2]==rpy[3]; rpy[3]==rpy[4]; rpy[4]==rpy[5]; rpy[5]==rpy[6]; rpy[6]==rpy[7]; rpy[8]==rpy[9]; rpy[9]==rpy[10]; rpy[10]==rpy[11]; }
		if (rpy[0]==' '){rpy[0]==rpy[1]; rpy[1]==rpy[2]; rpy[2]==rpy[3]; rpy[3]==rpy[4]; rpy[4]==rpy[5]; rpy[5]==rpy[6]; rpy[6]==rpy[7]; rpy[8]==rpy[9]; rpy[9]==rpy[10]; rpy[10]==rpy[11]; }
		if (rpy[0]==' '){rpy[0]==rpy[1]; rpy[1]==rpy[2]; rpy[2]==rpy[3]; rpy[3]==rpy[4]; rpy[4]==rpy[5]; rpy[5]==rpy[6]; rpy[6]==rpy[7]; rpy[8]==rpy[9]; rpy[9]==rpy[10]; rpy[10]==rpy[11]; }
		if (rpy[0]==' '){rpy[0]==rpy[1]; rpy[1]==rpy[2]; rpy[2]==rpy[3]; rpy[3]==rpy[4]; rpy[4]==rpy[5]; rpy[5]==rpy[6]; rpy[6]==rpy[7]; rpy[8]==rpy[9]; rpy[9]==rpy[10]; rpy[10]==rpy[11]; }
		if (rpz[0]==' '){rpz[0]==rpz[1]; rpz[1]==rpz[2]; rpz[2]==rpz[3]; rpz[3]==rpz[4]; rpz[4]==rpz[5]; rpz[5]==rpz[6]; rpz[6]==rpz[7]; rpz[8]==rpz[9]; rpz[9]==rpz[10]; rpz[10]==rpz[11]; }
		if (rpz[0]==' '){rpz[0]==rpz[1]; rpz[1]==rpz[2]; rpz[2]==rpz[3]; rpz[3]==rpz[4]; rpz[4]==rpz[5]; rpz[5]==rpz[6]; rpz[6]==rpz[7]; rpz[8]==rpz[9]; rpz[9]==rpz[10]; rpz[10]==rpz[11]; }
		if (rpz[0]==' '){rpz[0]==rpz[1]; rpz[1]==rpz[2]; rpz[2]==rpz[3]; rpz[3]==rpz[4]; rpz[4]==rpz[5]; rpz[5]==rpz[6]; rpz[6]==rpz[7]; rpz[8]==rpz[9]; rpz[9]==rpz[10]; rpz[10]==rpz[11]; }
		if (rpz[0]==' '){rpz[0]==rpz[1]; rpz[1]==rpz[2]; rpz[2]==rpz[3]; rpz[3]==rpz[4]; rpz[4]==rpz[5]; rpz[5]==rpz[6]; rpz[6]==rpz[7]; rpz[8]==rpz[9]; rpz[9]==rpz[10]; rpz[10]==rpz[11]; }
		//elimina erros dos vetores da raiz do ligante
		if (rlx[0]==' '){rlx[0]==rlx[1]; rlx[1]==rlx[2]; rlx[2]==rlx[3]; rlx[3]==rlx[4]; rlx[4]==rlx[5]; rlx[5]==rlx[6]; rlx[6]==rlx[7]; rlx[8]==rlx[9]; rlx[9]==rlx[10]; rlx[10]==rlx[11]; }
		if (rlx[0]==' '){rlx[0]==rlx[1]; rlx[1]==rlx[2]; rlx[2]==rlx[3]; rlx[3]==rlx[4]; rlx[4]==rlx[5]; rlx[5]==rlx[6]; rlx[6]==rlx[7]; rlx[8]==rlx[9]; rlx[9]==rlx[10]; rlx[10]==rlx[11]; }
		if (rlx[0]==' '){rlx[0]==rlx[1]; rlx[1]==rlx[2]; rlx[2]==rlx[3]; rlx[3]==rlx[4]; rlx[4]==rlx[5]; rlx[5]==rlx[6]; rlx[6]==rlx[7]; rlx[8]==rlx[9]; rlx[9]==rlx[10]; rlx[10]==rlx[11]; }
		if (rlx[0]==' '){rlx[0]==rlx[1]; rlx[1]==rlx[2]; rlx[2]==rlx[3]; rlx[3]==rlx[4]; rlx[4]==rlx[5]; rlx[5]==rlx[6]; rlx[6]==rlx[7]; rlx[8]==rlx[9]; rlx[9]==rlx[10]; rlx[10]==rlx[11]; }
		if (rly[0]==' '){rly[0]==rly[1]; rly[1]==rly[2]; rly[2]==rly[3]; rly[3]==rly[4]; rly[4]==rly[5]; rly[5]==rly[6]; rly[6]==rly[7]; rly[8]==rly[9]; rly[9]==rly[10]; rly[10]==rly[11]; }
		if (rly[0]==' '){rly[0]==rly[1]; rly[1]==rly[2]; rly[2]==rly[3]; rly[3]==rly[4]; rly[4]==rly[5]; rly[5]==rly[6]; rly[6]==rly[7]; rly[8]==rly[9]; rly[9]==rly[10]; rly[10]==rly[11]; }
		if (rly[0]==' '){rly[0]==rly[1]; rly[1]==rly[2]; rly[2]==rly[3]; rly[3]==rly[4]; rly[4]==rly[5]; rly[5]==rly[6]; rly[6]==rly[7]; rly[8]==rly[9]; rly[9]==rly[10]; rly[10]==rly[11]; }
		if (rly[0]==' '){rly[0]==rly[1]; rly[1]==rly[2]; rly[2]==rly[3]; rly[3]==rly[4]; rly[4]==rly[5]; rly[5]==rly[6]; rly[6]==rly[7]; rly[8]==rly[9]; rly[9]==rly[10]; rly[10]==rly[11]; }
		if (rlz[0]==' '){rlz[0]==rlz[1]; rlz[1]==rlz[2]; rlz[2]==rlz[3]; rlz[3]==rlz[4]; rlz[4]==rlz[5]; rlz[5]==rlz[6]; rlz[6]==rlz[7]; rlz[8]==rlz[9]; rlz[9]==rlz[10]; rlz[10]==rlz[11]; }
		if (rlz[0]==' '){rlz[0]==rlz[1]; rlz[1]==rlz[2]; rlz[2]==rlz[3]; rlz[3]==rlz[4]; rlz[4]==rlz[5]; rlz[5]==rlz[6]; rlz[6]==rlz[7]; rlz[8]==rlz[9]; rlz[9]==rlz[10]; rlz[10]==rlz[11]; }
		if (rlz[0]==' '){rlz[0]==rlz[1]; rlz[1]==rlz[2]; rlz[2]==rlz[3]; rlz[3]==rlz[4]; rlz[4]==rlz[5]; rlz[5]==rlz[6]; rlz[6]==rlz[7]; rlz[8]==rlz[9]; rlz[9]==rlz[10]; rlz[10]==rlz[11]; }
		if (rlz[0]==' '){rlz[0]==rlz[1]; rlz[1]==rlz[2]; rlz[2]==rlz[3]; rlz[3]==rlz[4]; rlz[4]==rlz[5]; rlz[5]==rlz[6]; rlz[6]==rlz[7]; rlz[8]==rlz[9]; rlz[9]==rlz[10]; rlz[10]==rlz[11]; }
		//elimina erros dos vetores do aceitador/doador da proteina
		if (bpx[0]==' '){bpx[0]==bpx[1]; bpx[1]==bpx[2]; bpx[2]==bpx[3]; bpx[3]==bpx[4]; bpx[4]==bpx[5]; bpx[5]==bpx[6]; bpx[6]==bpx[7]; bpx[8]==bpx[9]; bpx[9]==bpx[10]; bpx[10]==bpx[11]; }
		if (bpx[0]==' '){bpx[0]==bpx[1]; bpx[1]==bpx[2]; bpx[2]==bpx[3]; bpx[3]==bpx[4]; bpx[4]==bpx[5]; bpx[5]==bpx[6]; bpx[6]==bpx[7]; bpx[8]==bpx[9]; bpx[9]==bpx[10]; bpx[10]==bpx[11]; }
		if (bpx[0]==' '){bpx[0]==bpx[1]; bpx[1]==bpx[2]; bpx[2]==bpx[3]; bpx[3]==bpx[4]; bpx[4]==bpx[5]; bpx[5]==bpx[6]; bpx[6]==bpx[7]; bpx[8]==bpx[9]; bpx[9]==bpx[10]; bpx[10]==bpx[11]; }
		if (bpx[0]==' '){bpx[0]==bpx[1]; bpx[1]==bpx[2]; bpx[2]==bpx[3]; bpx[3]==bpx[4]; bpx[4]==bpx[5]; bpx[5]==bpx[6]; bpx[6]==bpx[7]; bpx[8]==bpx[9]; bpx[9]==bpx[10]; bpx[10]==bpx[11]; }
		if (bpy[0]==' '){bpy[0]==bpy[1]; bpy[1]==bpy[2]; bpy[2]==bpy[3]; bpy[3]==bpy[4]; bpy[4]==bpy[5]; bpy[5]==bpy[6]; bpy[6]==bpy[7]; bpy[8]==bpy[9]; bpy[9]==bpy[10]; bpy[10]==bpy[11]; }
		if (bpy[0]==' '){bpy[0]==bpy[1]; bpy[1]==bpy[2]; bpy[2]==bpy[3]; bpy[3]==bpy[4]; bpy[4]==bpy[5]; bpy[5]==bpy[6]; bpy[6]==bpy[7]; bpy[8]==bpy[9]; bpy[9]==bpy[10]; bpy[10]==bpy[11]; }
		if (bpy[0]==' '){bpy[0]==bpy[1]; bpy[1]==bpy[2]; bpy[2]==bpy[3]; bpy[3]==bpy[4]; bpy[4]==bpy[5]; bpy[5]==bpy[6]; bpy[6]==bpy[7]; bpy[8]==bpy[9]; bpy[9]==bpy[10]; bpy[10]==bpy[11]; }
		if (bpy[0]==' '){bpy[0]==bpy[1]; bpy[1]==bpy[2]; bpy[2]==bpy[3]; bpy[3]==bpy[4]; bpy[4]==bpy[5]; bpy[5]==bpy[6]; bpy[6]==bpy[7]; bpy[8]==bpy[9]; bpy[9]==bpy[10]; bpy[10]==bpy[11]; }
		if (bpz[0]==' '){bpz[0]==bpz[1]; bpz[1]==bpz[2]; bpz[2]==bpz[3]; bpz[3]==bpz[4]; bpz[4]==bpz[5]; bpz[5]==bpz[6]; bpz[6]==bpz[7]; bpz[8]==bpz[9]; bpz[9]==bpz[10]; bpz[10]==bpz[11]; }
		if (bpz[0]==' '){bpz[0]==bpz[1]; bpz[1]==bpz[2]; bpz[2]==bpz[3]; bpz[3]==bpz[4]; bpz[4]==bpz[5]; bpz[5]==bpz[6]; bpz[6]==bpz[7]; bpz[8]==bpz[9]; bpz[9]==bpz[10]; bpz[10]==bpz[11]; }
		if (bpz[0]==' '){bpz[0]==bpz[1]; bpz[1]==bpz[2]; bpz[2]==bpz[3]; bpz[3]==bpz[4]; bpz[4]==bpz[5]; bpz[5]==bpz[6]; bpz[6]==bpz[7]; bpz[8]==bpz[9]; bpz[9]==bpz[10]; bpz[10]==bpz[11]; }
		if (bpz[0]==' '){bpz[0]==bpz[1]; bpz[1]==bpz[2]; bpz[2]==bpz[3]; bpz[3]==bpz[4]; bpz[4]==bpz[5]; bpz[5]==bpz[6]; bpz[6]==bpz[7]; bpz[8]==bpz[9]; bpz[9]==bpz[10]; bpz[10]==bpz[11]; }
		//elimina erros dos vetores do aceitador/doador do ligante
		if (blx[0]==' '){blx[0]==blx[1]; blx[1]==blx[2]; blx[2]==blx[3]; blx[3]==blx[4]; blx[4]==blx[5]; blx[5]==blx[6]; blx[6]==blx[7]; blx[8]==blx[9]; blx[9]==blx[10]; blx[10]==blx[11]; }
		if (blx[0]==' '){blx[0]==blx[1]; blx[1]==blx[2]; blx[2]==blx[3]; blx[3]==blx[4]; blx[4]==blx[5]; blx[5]==blx[6]; blx[6]==blx[7]; blx[8]==blx[9]; blx[9]==blx[10]; blx[10]==blx[11]; }
		if (blx[0]==' '){blx[0]==blx[1]; blx[1]==blx[2]; blx[2]==blx[3]; blx[3]==blx[4]; blx[4]==blx[5]; blx[5]==blx[6]; blx[6]==blx[7]; blx[8]==blx[9]; blx[9]==blx[10]; blx[10]==blx[11]; }
		if (blx[0]==' '){blx[0]==blx[1]; blx[1]==blx[2]; blx[2]==blx[3]; blx[3]==blx[4]; blx[4]==blx[5]; blx[5]==blx[6]; blx[6]==blx[7]; blx[8]==blx[9]; blx[9]==blx[10]; blx[10]==blx[11]; }
		if (bly[0]==' '){bly[0]==bly[1]; bly[1]==bly[2]; bly[2]==bly[3]; bly[3]==bly[4]; bly[4]==bly[5]; bly[5]==bly[6]; bly[6]==bly[7]; bly[8]==bly[9]; bly[9]==bly[10]; bly[10]==bly[11]; }
		if (bly[0]==' '){bly[0]==bly[1]; bly[1]==bly[2]; bly[2]==bly[3]; bly[3]==bly[4]; bly[4]==bly[5]; bly[5]==bly[6]; bly[6]==bly[7]; bly[8]==bly[9]; bly[9]==bly[10]; bly[10]==bly[11]; }
		if (bly[0]==' '){bly[0]==bly[1]; bly[1]==bly[2]; bly[2]==bly[3]; bly[3]==bly[4]; bly[4]==bly[5]; bly[5]==bly[6]; bly[6]==bly[7]; bly[8]==bly[9]; bly[9]==bly[10]; bly[10]==bly[11]; }
		if (bly[0]==' '){bly[0]==bly[1]; bly[1]==bly[2]; bly[2]==bly[3]; bly[3]==bly[4]; bly[4]==bly[5]; bly[5]==bly[6]; bly[6]==bly[7]; bly[8]==bly[9]; bly[9]==bly[10]; bly[10]==bly[11]; }
		if (blz[0]==' '){blz[0]==blz[1]; blz[1]==blz[2]; blz[2]==blz[3]; blz[3]==blz[4]; blz[4]==blz[5]; blz[5]==blz[6]; blz[6]==blz[7]; blz[8]==blz[9]; blz[9]==blz[10]; blz[10]==blz[11]; }
		if (blz[0]==' '){blz[0]==blz[1]; blz[1]==blz[2]; blz[2]==blz[3]; blz[3]==blz[4]; blz[4]==blz[5]; blz[5]==blz[6]; blz[6]==blz[7]; blz[8]==blz[9]; blz[9]==blz[10]; blz[10]==blz[11]; }
		if (blz[0]==' '){blz[0]==blz[1]; blz[1]==blz[2]; blz[2]==blz[3]; blz[3]==blz[4]; blz[4]==blz[5]; blz[5]==blz[6]; blz[6]==blz[7]; blz[8]==blz[9]; blz[9]==blz[10]; blz[10]==blz[11]; }
		if (blz[0]==' '){blz[0]==blz[1]; blz[1]==blz[2]; blz[2]==blz[3]; blz[3]==blz[4]; blz[4]==blz[5]; blz[5]==blz[6]; blz[6]==blz[7]; blz[8]==blz[9]; blz[9]==blz[10]; blz[10]==blz[11]; }
		// transforma os 4 vetores em numeros
		// proteina
		//RP = C1
		num_rpx = atof(rpx);
		num_rpy = atof(rpy);		 		 
		num_rpz = atof(rpz);
		//ADP = A
		num_bpx = atof(bpx);
		num_bpy = atof(bpy);
		num_bpz = atof(bpz);
		// ligante
		//RL = C2
		num_rlx = atof(rlx);
		num_rly = atof(rly);
		num_rlz = atof(rlz);
		//ADL = B
		num_blx = atof(blx);
		num_bly = atof(bly);
		num_blz = atof(blz);
		//--------------------------
		//cout << num_blx << endl;						
		//-----------------A-B----X
		double ABXt1;
		ABXt1 = (num_blx-num_bpx);
		double ABXt2;
		ABXt2 = (ABXt1*ABXt1);
		//-------------A-B-------Y
		double ABYt1;
		ABYt1 = (num_bly-num_bpy);
		double ABYt2;
		ABYt2 = (ABYt1*ABYt1);
		//-------------A-B-------Z
		double ABZt1;
		ABZt1 = (num_blz-num_bpz);
		double ABZt2;
		ABZt2 = (ABZt1*ABZt1);
		//-------soma AB
		double ABt1;
		ABt1 = ABXt2+ABYt2;
		double ABt2;
		ABt2 = ABt1+ABZt2;
		//--------D entre a e b
		double DAB;
		DAB = sqrt(ABt2);
		//----------------------------------ANGULO-UM-------------------------
		//-------------A-C1-----------X
		double AC1Xt1;
		AC1Xt1 = (num_rpx-num_bpx);
		double AC1Xt2;
		AC1Xt2 = (AC1Xt1*AC1Xt1);
		//-------------A-C1-------Y
		double AC1Yt1;
		AC1Yt1 = (num_rpy-num_bpy);
		double AC1Yt2;
		AC1Yt2 = (AC1Yt1*AC1Yt1);
		//-------------A-C1-------Z
		double AC1Zt1;
		AC1Zt1 = (num_rpz-num_bpz);
		double AC1Zt2;
		AC1Zt2 = (AC1Zt1*AC1Zt1);
		//-------soma AC1
		double AC1t1;
		AC1t1 = AC1Xt2+AC1Yt2;
		double AC1t2;
		AC1t2 = AC1t1+AC1Zt2;
		//--------D entre a e C1
		double DAC1;
		DAC1 = sqrt(AC1t2);
		//________________________________________________________
		//-------------B-C1-----------X
		double BC1Xt1;
		BC1Xt1 = (num_rpx-num_blx);
		double BC1Xt2;
		BC1Xt2 = (BC1Xt1*BC1Xt1);
		//-------------B-C1-------Y
		double BC1Yt1;
		BC1Yt1 = (num_rpy-num_bly);
		double BC1Yt2;
		BC1Yt2 = (BC1Yt1*BC1Yt1);
		//-------------B-C1-------Z
		double BC1Zt1;
		BC1Zt1 = (num_rpz-num_blz);
		double BC1Zt2;
		BC1Zt2 = (BC1Zt1*BC1Zt1);
		//-------soma BC1
		double BC1t1;
		BC1t1 = BC1Xt2+BC1Yt2;
		double BC1t2;
		BC1t2 = BC1t1+BC1Zt2;
		//--------D entre a e C1
		double DBC1;
		DBC1 = sqrt(BC1t2);
		//_______________________________________________
		//----------------------------------ANGULO-DOIS-------------------------
		//-------------A-C2-----------X
		double AC2Xt1;
		AC2Xt1 = (num_rlx-num_bpx);
		double AC2Xt2;
		AC2Xt2 = (AC2Xt1*AC2Xt1);
		//-------------A-C2-------Y
		double AC2Yt1;
		AC2Yt1 = (num_rly-num_bpy);
		double AC2Yt2;
		AC2Yt2 = (AC2Yt1*AC2Yt1);
		//-------------A-C2-------Z
		double AC2Zt1;
		AC2Zt1 = (num_rlz-num_bpz);
		double AC2Zt2;
		AC2Zt2 = (AC2Zt1*AC2Zt1);
		//-------soma AC2
		double AC2t1;
		AC2t1 = AC2Xt2+AC2Yt2;
		double AC2t2;
		AC2t2 = AC2t1+AC2Zt2;
		//--------D entre a e C2
		double DAC2;
		DAC2 = sqrt(AC2t2);
		//________________________________________________________
		//-------------B-C2-----------X
		double BC2Xt1;
		BC2Xt1 = (num_rlx-num_blx);
		double BC2Xt2;
		BC2Xt2 = (BC2Xt1*BC2Xt1);
		//-------------B-C2-------Y
		double BC2Yt1;
		BC2Yt1 = (num_rly-num_bly);
		double BC2Yt2;
		BC2Yt2 = (BC2Yt1*BC2Yt1);
		//-------------B-C2-------Z
		double BC2Zt1;
		BC2Zt1 = (num_rlz-num_blz);
		double BC2Zt2;
		BC2Zt2 = (BC2Zt1*BC2Zt1);
		//-------soma BC2
		double BC2t1;
		BC2t1 = BC2Xt2+BC2Yt2;
		double BC2t2;
		BC2t2 = BC2t1+BC2Zt2;
		//--------D entre a e C2
		double DBC2;
		DBC2 = sqrt(BC2t2);
		//_______________________________________________
		//----------------ANGULOS
		double DAB_2; //DAB ao quadrado
		DAB_2 = (DAB*DAB);
		//-----------ANGULO-----UM
		double DAC1_2;
		DAC1_2 = (DAC1*DAC1);
		double DBC1_2;
		DBC1_2 = (DBC1*DBC1);
		//soma AB AC1
		double S_DAB_DAC1_DBC1_t1; //DAB + DAC1
		S_DAB_DAC1_DBC1_t1 = (DAB_2+DAC1_2);
		double S_DAB_DAC1_DBC1_t2; //DAB + DAC1 - DBC1
		S_DAB_DAC1_DBC1_t2 = (S_DAB_DAC1_DBC1_t1-DBC1_2);
		//AB * AC1
		double r2V_DAB_DAC1_t1;
		r2V_DAB_DAC1_t1 = (DAB*DAC1);
		double r2V_DAB_DAC1_t2;
		r2V_DAB_DAC1_t2 = (2*r2V_DAB_DAC1_t1);
		//--------final------angulo 1
		double r2V_DAB_DAC1_t3;
		r2V_DAB_DAC1_t3 = (S_DAB_DAC1_DBC1_t2/r2V_DAB_DAC1_t2);
		double angulo_1_t1;
		angulo_1_t1 = acos(r2V_DAB_DAC1_t3);
		//-----radianos em graus
		double angulo_1_t2;
		angulo_1_t2 = (angulo_1_t1*180);
		double angulo_1;
		angulo_1 = (angulo_1_t2/3.141618);
		//----------------------------------
		//----------------------------------ANGULO-DOIS
		//-----------ANGULO-----Dois
		double DAC2_2;
		DAC2_2 = (DAC2*DAC2);
		double DBC2_2;
		DBC2_2 = (DBC2*DBC2);
		//soma AB AC1
		double S_DAB_DAC2_DBC2_t1; //DAB + DAC2
		S_DAB_DAC2_DBC2_t1 = (DAB_2+DBC2_2);
		double S_DAB_DAC2_DBC2_t2; //DAB + DAC1 - DBC2
		S_DAB_DAC2_DBC2_t2 = (S_DAB_DAC2_DBC2_t1-DAC2_2);
		//AB * AC2
		double r2V_DAB_DAC2_t1;
		r2V_DAB_DAC2_t1 = (DAB*DBC2);
		double r2V_DAB_DAC2_t2;
		r2V_DAB_DAC2_t2 = (2*r2V_DAB_DAC2_t1);
		//--------final------angulo 2
		double r2V_DAB_DAC2_t3;
		r2V_DAB_DAC2_t3 = (S_DAB_DAC2_DBC2_t2/r2V_DAB_DAC2_t2);
		double angulo_2_t1;
		angulo_2_t1 = acos(r2V_DAB_DAC2_t3);
		//-----radianos em graus
		double angulo_2_t2;
		angulo_2_t2 = (angulo_2_t1*180);
		double angulo_2;
		angulo_2 = (angulo_2_t2/3.141618);
		//----------------------------------
		//cout << "angulo 1   = " << angulo_1 << endl;
		//cout << "angulo 2   = " << angulo_2 << endl;
		//cout << endl;
		//--------------- SALVA ARQUIVO PDB COM ATOMOS E RAIZES OS ATOMOS MAIS PROXIMOS
		//int num1_rl_atomo = 0;
		int num1_bl_atomo = 0;
		//int num1_rp_atomo = 0;
		int num1_bp_atomo = 0;
		//int num_rl_atomo = 0;
		int num_bl_atomo = 0;
		//int num_rp_atomo = 0;
		int num_bp_atomo = 0;
		//num_rl_atomo = atoi(rl_atomo);
		num_bl_atomo = atoi(bl_atomo);
		//num_rp_atomo = atoi(rp_atomo);
		num_bp_atomo = atoi(bp_atomo);
		if ((angulo_1>=60) && (angulo_2>=60))
		{
			//contador de ligacoes de hidrogenio
			bond_count++;
			//--------------saida-ligante-TEMPORARIA
			//if (num_bl_atomo>num1_bl_atomo)
			//{
			num1_bl_atomo = num_bl_atomo;
			l_result << "HETATM" << bl_atomo[2] << bl_atomo[3] << bl_atomo[4] << bl_atomo[5] << bl_atomo[6] << bl_atomo[7] << " " << bl_nome_atomo << "BLK" << bl_resto << "     " << c1_blx[0] << c1_blx[1] << c1_blx[2] << c1_blx[3] << c1_blx[4] << c1_blx[5] << c1_blx[6] << c1_blx[7] << c1_bly[0] << c1_bly[1] << c1_bly[2] << c1_bly[3] << c1_bly[4] << c1_bly[5] << c1_bly[6] << c1_bly[7] << c1_blz[0] << c1_blz[1] << c1_blz[2] << c1_blz[3] << c1_blz[4] << c1_blz[5] << c1_blz[6] << c1_blz[7] <<  "                         " << bl_enter;
			//}
			//------------------------------------------
			//--------------saida-proteina-TEMPORARIA
			//if (num_bp_atomo>num1_bp_atomo)
			//{
			num1_bp_atomo = num_bp_atomo;
			p_result << "ATOM " << bp_atomo  << " " << bp_nome_atomo << bp_aminoacido << bp_num_aminoacido << "   " << c1_bpx << c1_bpy << c1_bpz << "                         " << bp_enter;
			//}
			//------------------------------------------
			char distance[14] = "             ";
			sprintf(distance,"%.3f",DAB);
			strcat(distance, " ");
			int k=0;
			for (k = 0; k < 11; k++)
			{if (distance[k]!='\0'){distances << distance[k];}}
			distances << endl;
		}//-----------if angulo1 e 2 obedecem os limites
	} 	 
	p_result.close(); //-----fecha saida proteina
	l_result.close(); //-----fecha saida ligante
	proteina.close();
	ligante.close(); 	 
	distances.close();
	// mostra a contagem total de pontes de hidrogenio
	/*
	cout << endl << bond_count << endl;
	ofstream bound_total;
	bound_total.open("bond_count.txt");
	bound_total << bond_count;
	bound_total.close();
	*/
	//------------------saida: resultados detalhados
	//entradas da proteina e ligante
	char p_line[82];
	char l_line[82];
	char b_distances[12];
	//gera arquivo de log
	ofstream bonds_log;
	bonds_log.open("bonds.log");
	//abre arquivo para depois comparar o limite maximo de ligacoes de hidrogenio
	ofstream limit_l;
	ofstream limit_p;
	//para proteina
	limit_p.open("limit_p.tmp");
	//para ligante
	limit_l.open("limit_l.tmp");
	bonds_log << "Hydrogen bonds found: " << bond_count << endl << endl;
	//cout << "Hydrogen bonds found: " << bond_count << endl << endl; 
	hydrogen_B = bond_count;
	ifstream p_bonds;
	ifstream l_bonds;
	p_bonds.open("p_result.tmp");
	l_bonds.open("l_result.tmp");
	//d_bonds.open("dist_result.tmp");
	while (!p_bonds.eof() || !l_bonds.eof())
	{
		p_bonds.get(p_line, sizeof(p_line), '\0');
		l_bonds.get(l_line, sizeof(l_line), '\0');
		//d_bonds.get(b_distances, sizeof(b_distances), '\0');
		bonds_log << p_line << l_line;
		// << "distance: " << b_distances;
		bonds_log << endl;
	} 
	p_bonds.close();
	l_bonds.close();
	bonds_log << "--------------------------------------------------------" << endl;
	bonds_log << "        TABLE: INTERMOLECULAR HYDROGEN BONDS            " << endl;
	bonds_log << endl;
	bonds_log << "       Protein               Ligand          Distance(A)" << endl;
	ifstream d_bonds;
	ifstream p_bonds_final;
	ifstream l_bonds_final;
	p_bonds_final.open("p_result.tmp");
	l_bonds_final.open("l_result.tmp");
	d_bonds.open("dist_result.tmp");
	while (!p_bonds_final.eof() || !l_bonds_final.eof() || !d_bonds.eof())
	{
		p_bonds_final.get(p_line, sizeof(p_line), '\0');
		l_bonds_final.get(l_line, sizeof(l_line), '\0');
		d_bonds.get(b_distances, sizeof(b_distances), '\0');
		//proteina nome aminoacido
		bonds_log << p_line[16] << p_line[17] << p_line[18] << p_line[19] << p_line[20];
		limit_p << p_line[16] << p_line[17] << p_line[18] << p_line[19] << p_line[20];
		//proteina numero aminoacido
		bonds_log << p_line[21] << p_line[22] << p_line[23] << p_line[24] << p_line[25] << p_line[26];
		limit_p << p_line[21] << p_line[22] << p_line[23] << p_line[24] << p_line[25] << p_line[26];
		//proteina nome atomo
		bonds_log << p_line[12] << p_line[13] << p_line[14] << p_line[15];
		limit_p << p_line[12] << p_line[13] << p_line[14] << p_line[15];
		//proteina numero atomo
		bonds_log << p_line[6] << p_line[7] << p_line[8] << p_line[9] << p_line[10] << p_line[11];
		limit_p << p_line[6] << p_line[7] << p_line[8] << p_line[9] << p_line[10] << p_line[11] << endl;
		//proteina divisoria
		bonds_log << "   ";
		//----------------------------------------------
		//ligante nome aminoacido
		bonds_log << l_line[16] << l_line[17] << l_line[18] << l_line[19] << l_line[20];
		limit_l << l_line[16] << l_line[17] << l_line[18] << l_line[19] << l_line[20];
		//ligante numero aminoacido
		//bonds_log << l_line[22] << l_line[23] << l_line[24] << l_line[25] << l_line[26];
		//ligante nome atomo
		bonds_log << l_line[12] << l_line[13] << l_line[14] << l_line[15];
		limit_l << l_line[12] << l_line[13] << l_line[14] << l_line[15];
		//ligante numero atomo
		bonds_log << l_line[6] << l_line[7] << l_line[8] << l_line[9] << l_line[10] << l_line[11];
		limit_l << l_line[6] << l_line[7] << l_line[8] << l_line[9] << l_line[10] << l_line[11] << endl;
		//liagante  divisoria
		bonds_log << "         ";
		bonds_log << b_distances;
		//bonds_log << endl;
	} 
	p_bonds_final.close();
	l_bonds_final.close();
	d_bonds.close();
	limit_p.close();
	limit_l.close();
	bonds_log.close();
}
//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
//-----------------------------------------------------------------
//---------------SAIDA-FORMATO-PDB------------------
//--------------------------------------------------------------
void saida_PDB()
{
	//saida final
	vector<string> v;
	v.clear();
	split(ligand_name, '.', v);
	char hb_file_name[100];
	strcpy(hb_file_name, v[0].c_str());
	strcat(hb_file_name, "_H-Bonds.pdb");
	ofstream result_PDB;
	result_PDB.open(hb_file_name); 
	//entrada (TEMP) 	 
	ifstream p_result;
	p_result.open("p_result.tmp");
	ifstream l_result;
	l_result.open("l_result.tmp");
	char ATOM[7];
	char esp1[2];
	char esp2[4];
	char esp3[26];
	char bp_enter[3];
	//aceitador/doador
	char bp_atomo[7];
	char bp_nome_atomo[5];
	char bp_aminoacido[6];
	char bp_num_aminoacido[6];
	int num_bp_num_aminoacido = 0;
	int num_bp_atomo = 0;
	int atomo_anterior = 0;
	int num_p_num_aminoacido = 0;
	int num_p_atomo = 0;
	char bpx[9];
	char bpy[9];
	char bpz[9];
	//pega o aminoacido pertencente aos atomos envolvidos nas ligacoes de H
	//variaveis do ifstream protein
	char p_atomo[7];
	char p_nome_atomo[5];
	char p_aminoacido[6];
	char p_num_aminoacido[6];
	char px[9];
	char py[9];
	char pz[9];
	char p_enter[3];
	int num_aminoacido_anterior = 0;
	while (!p_result.eof())
	{
		p_result.get(ATOM,sizeof(ATOM),'\0');
		p_result.get(bp_atomo,sizeof(bp_atomo),'\0');
		p_result.get(esp1,sizeof(esp1),'\0');
		p_result.get(bp_nome_atomo,sizeof(bp_nome_atomo),'\0');
		p_result.get(bp_aminoacido,sizeof(bp_aminoacido),'\0');
		p_result.get(bp_num_aminoacido,sizeof(bp_num_aminoacido),'\0');
		p_result.get(esp2,sizeof(esp2),'\0');
		p_result.get(bpx,sizeof(bpx),'\0');
		p_result.get(bpy,sizeof(bpy),'\0');
		p_result.get(bpz,sizeof(bpz),'\0');
		p_result.get(esp3,sizeof(esp3),'\0');    
		p_result.get(bp_enter,sizeof(bp_enter),'\0');  				 
		//cout << "*" << ATOM << "*" << bp_atomo << "*" << esp1 << "*" << bp_nome_atomo << "*" << bp_aminoacido << "*" << bp_num_aminoacido << "*" << esp2 << "*" << bpx << "*" << bpy << "*" << bpz << "*" << esp3 << "*" << bp_enter;
		num_bp_num_aminoacido = atoi(bp_num_aminoacido);
		num_bp_atomo = atoi(bp_atomo);
		ifstream protein;
		protein.open("protein.tmp");
		while (!protein.eof())
		{
			protein.get(p_atomo,sizeof(p_atomo),'\0');
			protein.ignore(1);
			protein.get(p_nome_atomo,sizeof(p_nome_atomo),'\0');
			protein.get(p_aminoacido,sizeof(p_aminoacido),'\0');
			protein.get(p_num_aminoacido,sizeof(p_num_aminoacido),'\0');
			protein.get(px,sizeof(px),'\0');
			protein.get(py,sizeof(py),'\0');
			protein.get(pz,sizeof(pz),'\0');
			protein.get(p_enter,sizeof(p_enter),'\0');
			//cout << "*" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << px << "*" << py << "*" << pz << "*" << p_enter;
			num_p_num_aminoacido = atoi(p_num_aminoacido);
			num_p_atomo = atoi(p_atomo);
			//cout << num_p_num_aminoacido << endl;
			//cout << num_p_atomo << endl;
			//cout << p_num_aminoacido << endl;
			if (num_p_num_aminoacido==num_bp_num_aminoacido/* && num_aminoacido_anterior!=num_p_num_aminoacido*/)
			{
				if(num_p_atomo<num_bp_atomo){result_PDB << "ATOM "  << p_atomo << " " << esp1 << p_nome_atomo << p_aminoacido << p_num_aminoacido << esp2 << px << py << pz  << esp3 << p_enter;
				//cout << "*" << "ATOM " << "*" << p_atomo << " " << "*" << esp1 << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << esp2 << "*" << px << "*" << py << "*" << pz << "*" << esp3 << "*" << p_enter;
				}
				if(num_p_atomo==num_bp_atomo){result_PDB << ATOM << bp_atomo << esp1 << bp_nome_atomo << bp_aminoacido << bp_num_aminoacido << esp2 << bpx << bpy << bpz << esp3 << bp_enter;
				//cout << "*" << ATOM << "*" << bp_atomo << "*" << esp1 << "*" << bp_nome_atomo << "*" << bp_aminoacido << "*" << bp_num_aminoacido << "*" << esp2 << "*" << bpx << "*" << bpy << "*" << bpz << "*" << esp3 << "*" << bp_enter;
				}
				if(num_p_atomo>num_bp_atomo){result_PDB<< "ATOM " << p_atomo << " " << esp1 << p_nome_atomo << p_aminoacido << p_num_aminoacido << esp2 << px << py << pz << esp3 << p_enter;
				//cout << "*" << "ATOM " << "*" << p_atomo << " " << "*" << esp1 << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << esp2 << "*" << px << "*" << py << "*" << pz << "*" << esp3 << "*" << p_enter;
				}
			}//if
		}//while protein
		protein.close();
		num_aminoacido_anterior = num_bp_num_aminoacido;
	}//while p_result
	char ch[2];
	while (!l_result.eof())
	{
		l_result.get(ch, sizeof(ch),'\0');
		result_PDB << ch;
	}//while l_result
	result_PDB.close();
	l_result.close();
	p_result.close();
}//END
//-------------------------------------------------
//FUNCAO QUE SALVA ARQUIVO DO LIGANTE INTEIRO NO PDB
//-------------------------------------------------
void salva_ligante()
{
	char num_atomo[7];
	char nome_atomo[5];
	char x[9];
	char y[9];
	char z[9];
	char resto1[7];
	ifstream ligand;
	ligand.open("ligand.tmp");
	ofstream l_result;
	l_result.open("l_result.tmp");
	while (!ligand.eof())
	{
		ligand.ignore(2); 	   
		ligand.get(num_atomo, sizeof(num_atomo), '\0');
		ligand.get(nome_atomo, sizeof(nome_atomo), '\0');
		ligand.ignore(5);
		ligand.get(x, sizeof(x), '\0');
		ligand.ignore(2);
		ligand.get(y, sizeof(y), '\0');
		ligand.ignore(2);
		ligand.get(z, sizeof(z), '\0');
		ligand.ignore(2);
		ligand.get(resto1, sizeof(resto1), '\0');	
		if (!ligand.eof())
		{		
			l_result << "HETATM" << num_atomo << " " << nome_atomo << "BLK" << "          " << x << y << z << "                          " << "\n";		
		}
		if (ligand.eof())
		{		
			l_result << "END" << "                                                                             " << "\n";
		}
	}			 
	ligand.close();
	l_result.close();
}
//retira os aminoacidos repetidos da proteina
void del_rep(char (*vec)[TAM])
{
//ja achou o total de linhas e agora vai copiar todo arquivo para a string1
char aux[lines_total][82];
int i = 0;
int j = 0;
      for(i = 0; i < lines_total; i++)
            strcpy(aux[i],vec[i]);  /* copiamos a matriz a aux */
      for(i = 0; i < lines_total; i++)  
       { for(j = 0; j < lines_total; j++)   
           {if(!strcmp(vec[i],aux[j]) && i != j )  /*i != j */
                {strcpy(aux[i],"REP");}
           }    
       }               
lines_total2 = 0; 
lines_total2--;
       for(i = 0, j = 0; i < lines_total; i++)      /*modificam o original*/
                if(strcmp(aux[i],"REP"))    /* se e rep nao copiamos */
                    {strcpy(vec[j++],aux[i]);lines_total2++;}
//return;   /* nao faz falta
}
//FUNCAO QUE DELETA ARQUIVOS TEMPORARIOS
int salva_proteina() /* ha probarlo */
{
	int L = 0;
	lines_total = 0;
	ifstream p_result;
	vector<string> v;
	v.clear();
	split(ligand_name, '.', v);
	char hb_file_name[100];
	strcpy(hb_file_name, v[0].c_str());
	strcat(hb_file_name, "_H-Bonds.pdb");
	p_result.open(hb_file_name);
	char line[82];
	//conta quantas linhas tem
	while (!p_result.eof())
	{
		p_result.get(line, sizeof(line), '\0');
		lines_total++; 
	}
	p_result.close();
	//cout << lines_total;
	int i;   
	char vec[lines_total][TAM]; 
	ifstream p_result2;
	p_result2.open(hb_file_name);
	while (!p_result2.eof() && L <=lines_total)
	{
		p_result2.get(line, sizeof(line), '\0');
		strcpy(vec[L], line); 
		//cout << line << endl << L << endl;
		L++;	  
	}
	p_result2.close();
	//p_result2.get(vec[L], sizeof(vec));
	//L = 0;
	//for(L = 0; L < lines_total; L++)
	// cout << vec[L];
	L = 0;
	del_rep(vec);
	ofstream p_result3;
	strcpy(hb_file_name, v[0].c_str());
	strcat(hb_file_name, "_H-Bonds.pdb");
	p_result3.open (hb_file_name); 
	for(i = 0; i <= lines_total2; i++)
			//cout << vec[i]; 
			p_result3 << vec[i];
			//printf("%s\n",vec[i]);
	p_result3.close();
}
//funcao que limita o numero 
//maximo de ligacoes de hidrogenio 
//que vao aparecer no arquivo de sa�a final
void limit()
{
	//arquivo contendo a lista de atomos que participam de ligacoes
	//para proteina
	ifstream p_list;
	p_list.open("limit_p.tmp");
	//contadores
	int line_p_count = 0;
	int line_l_count = 0;
	int line_d_count = 0;
	int p_position = 0;
	int l_position = 0;
	int d_position = 0;
	int d_smaller = 0;
	int copied = 0;
	int compared = 0;
	//strings
	char atom_l;
	char atom_l2;
	char atom_p;
	char atom_p2;
	char distance;
	char distance2;
	double distance_num;
	double distance2_num;
	//maximo numero de ligacoes padrao
	int MAX_O = 4;
	int MAX_C = 4;
	int MAX_P = 4;
	int MAX_S = 4;
	int MAX_F = 4;
	int MAX_N = 4;
	int MAX = 0;
	char ch[23];
	int TOTAL = 0;
	int p_TOTAL = 0;
	while (!p_list.eof())
	{ 
		p_list.get(ch, sizeof(ch), '\0');
		if (ch[21]=='\n')
		{
			p_TOTAL++;
		//cout << ch << "*" << p_TOTAL << "*";
		}
		//cout << ch;
	} 
	p_list.clear();              // forget we hit the end of file
	p_list.seekg(0, ios::beg);   // move to the start of the file
	int line_p = 0;
	char ch_p[p_TOTAL][23];
	while (!p_list.eof() || line_p!=p_TOTAL)
	{ 
		p_list.get(ch, sizeof(ch), '\0');
		strcpy(ch_p[line_p], ch);
		//cout << ch << "*" << ch_p[line_p] << "*";
		line_p++; 
	} 
	line_p = 0;
	line_p_count = 0; // a linha que esta sendo comparada com as outras
	//strcmp (szKey,szInput) != 0 //falso, diferente
	//strncmp(cs , ct ,n) 	 
	p_list.close(); 	 
	ifstream p_root;
	p_root.open("P_ROOT_0.tmp");
	int p_liged = 0;
	char p_liged_list[p_TOTAL][3];
	int p_file_line = 0;
	char p_root_line[49];
	//--------------------------------------------
	p_function:
	//--------------------------------------------		   
	while (!p_root.eof()&&((p_file_line-1)<line_p_count))
	{
		p_root.get(p_root_line, sizeof(p_root_line), '\0');
		//cout << p_root_line << "*" << endl;
		//cout << "*" << p_root_line[8] << "*" << endl;
		if (p_root_line[0]=='R' && p_root_line[8]!='H')
		{
			p_liged++;
		}
		if (p_root_line[0]=='B')
		{
			p_file_line++; /*cout << "@" << p_file_line << "@";*/} 
		if ((p_root_line[0]=='B')&&(((p_file_line-1)!=line_p_count)||(p_root_line[2]!=ch_p[line_p_count][15])|| (p_root_line[3]!=ch_p[line_p_count][16])||(p_root_line[4]!=ch_p[line_p_count][17])||(p_root_line[5]!=ch_p[line_p_count][18])||(p_root_line[6]!=ch_p[line_p_count][19])))
		{
			p_liged = 0;
		}
	}
	sprintf(p_liged_list[line_p_count], "%i",p_liged);
	//cout << "numero de atomos ligados ao" << ch_p[line_p_count] << " = " << p_liged << "*" << p_liged_list[line_p_count] << endl;
	p_liged = 0;
	line_p_count++;
	p_file_line = 0;
	p_root.clear();              // forget we hit the end of file
	p_root.seekg(0, ios::beg);   // move to the start of the file
	if (line_p_count<p_TOTAL){goto p_function;}
	p_root.close();
	//-------------- ABRE LIGANTE ---------------------
	ifstream l_list;
	l_list.open("limit_l.tmp");
	int l_TOTAL = 0;
	char ch2[17];
	while (!l_list.eof())
	{ 
		l_list.get(ch2, sizeof(ch2), '\0');
		//cout << ch2 << "*";
		if (ch2[15]=='\n')
		{
			l_TOTAL++;
			//cout << ch2 << "*" << l_TOTAL << "*";
		}
		//cout << ch;
	} 
	l_list.clear();              // forget we hit the end of file
	l_list.seekg(0, ios::beg);   // move to the start of the file
	int line_l = 0;
	char ch_l[l_TOTAL][17];
	while (!l_list.eof() || line_l!=l_TOTAL)
	{ 
		l_list.get(ch2, sizeof(ch2), '\0');
		strcpy(ch_l[line_l], ch2);
		//cout << ch2 << "*" << ch_l[line_p] << "*";
		line_l++; 
	}  	 
	line_l = 0;
	line_l_count = 0; // a linha que esta sendo comparada com as outras
	//strcmp (szKey,szInput) != 0 //falso, diferente
	//strncmp(cs , ct ,n) 	 
	l_list.close();
	ifstream l_root;
	l_root.open("L_ROOT_0.tmp");
	int l_file_line = 0;
	int l_liged = 0;
	char l_liged_list[l_TOTAL][3];
	char l_root_line[55];
	//--------------------------------------------
	l_function:
	//--------------------------------------------		   
	while (!l_root.eof()&&((l_file_line-1)<line_l_count))
	{
		l_root.get(l_root_line, sizeof(l_root_line), '\0');
		//cout << l_root_line << "*" << endl;
		if ((l_root_line[0]=='R') && (l_root_line[48]!='H') && (l_root_line[48]!='h'))
		{l_liged++;}
		if (l_root_line[0]=='B'){l_file_line++; /*cout << "@" << l_file_line << "@";*/} 
		if ((l_root_line[0]=='B')&&(((l_file_line-1)!=line_l_count)||(l_root_line[3]!=ch_l[line_l_count][9])||(l_root_line[4]!=ch_l[line_l_count][10])||(l_root_line[5]!=ch_l[line_l_count][11])||(l_root_line[6]!=ch_l[line_l_count][12])||(l_root_line[7]!=ch_l[line_l_count][13])))
		{l_liged = 0;}
	}
	sprintf(l_liged_list[line_l_count], "%i",l_liged);
	//cout << "numero de atomos ligados ao" << ch_l[line_l_count] << " = " << l_liged << "*" << l_liged_list[line_l_count] << endl;
	l_liged = 0;
	line_l_count++;
	l_file_line = 0;
	l_root.clear();              // forget we hit the end of file
	l_root.seekg(0, ios::beg);   // move to the start of the file
	if (line_l_count<l_TOTAL){goto l_function;}
	l_root.close();
	//-------------------------abre distancias-------------------
	char dist_line[p_TOTAL][12];
	char ch3[12];
	int line_d = 0;
	ifstream d_list;
	d_list.open("dist_result.tmp");
	while (!d_list.eof() && line_d<p_TOTAL)
	{
		d_list.get(ch3, sizeof(ch3),'\0');
		strcpy(dist_line[line_d], ch3);
		//cout << "*" << dist_line[line_d] << "*" << endl;
		line_d++; 	  
	}
	line_d = 0;
	line_d_count = 0;
	//------------------------------------------------------------
	//-----------------------------COME� O REFINAMENTO------------------
	//-----------------------------------------------------------------
	line_l_count = 0;
	line_p_count = 0;
	int liged = 0; //numero de atomos ja ligados
	int found = 0;
	int string_line_count = 0;
	ifstream limit_type;
	limit_type.open("limit_type.tmp");
	char type[8];
	char atom[8];
	refined = 0; 
	while (line_p_count<p_TOTAL)
	{
		liged = atoi(p_liged_list[line_p_count]);
		if (ch_p[line_p_count][12]=='C')
		{MAX = (MAX_C - liged); /*cout << endl << "C" <<endl;*/}
		if (ch_p[line_p_count][12]=='N')
		{MAX = (MAX_N - liged); /*cout << endl << "N" <<endl;*/}
		if (ch_p[line_p_count][12]=='P')
		{MAX = (MAX_P - liged); /*cout << endl << "P" <<endl;*/}
		if (ch_p[line_p_count][12]=='O')
		{MAX = (MAX_O - liged); /*cout << endl << "O" <<endl;*/}
		if (ch_p[line_p_count][12]=='S')
		{MAX = (MAX_S - liged); /*cout << endl << "S" <<endl;*/}
		while (string_line_count<p_TOTAL)
		{
			if ((ch_p[line_p_count][12]!='X')&&(ch_p[line_p_count][12]==ch_p[string_line_count][12])&& (ch_p[line_p_count][13]==ch_p[string_line_count][13])&& (ch_p[line_p_count][14]==ch_p[string_line_count][14])&& (ch_p[line_p_count][15]==ch_p[string_line_count][15])&& (ch_p[line_p_count][16]==ch_p[string_line_count][16])&& (ch_p[line_p_count][17]==ch_p[string_line_count][17])&& (ch_p[line_p_count][18]==ch_p[string_line_count][18])&& (ch_p[line_p_count][19]==ch_p[string_line_count][19]))
			{
			found++; //encontrou um atomo igual ao da string selecionada
			if (found > MAX)
				{
					ch_p[string_line_count][12]='X'; //marca este atomo para excluir da procura
					ch_l[string_line_count][9]='X';
					dist_line[string_line_count][0]='X';
				}
			}
			string_line_count++;
		}
		//cout << endl << "F" << found << endl;
		string_line_count = 0; //reseta as linhas que vao serem comparadas com a selecionada
		found = 0; //reseta nmero de encontrados
		line_p_count++; //seleciona a proxima string a ser comparada com as outras
		//d_list.close();
	}
	line_p_count = 0;
	line_l_count = 0;
	//-------------------------------------------------------------
	//------refinamento ligante-------------------------------------------------------
	//---------------------------------------------------------------
	MAX = 0;
	while (line_l_count<l_TOTAL)
	{
		liged = atoi(l_liged_list[line_l_count]);
		if (ch_l[line_l_count][6]=='C')
		{MAX = (MAX_C-liged); /*cout << endl << "C" <<endl;*/}
		if (ch_l[line_l_count][6]=='N')
		{MAX = (MAX_N-liged); /*cout << endl << "N" <<endl;*/}
		if (ch_l[line_l_count][6]=='P')
		{MAX = (MAX_P-liged); /*cout << endl << "P" <<endl;*/}
		if (ch_l[line_l_count][6]=='O')
		{MAX = (MAX_O-liged); /*cout << endl << "O" <<endl;*/}
		if (ch_l[line_l_count][6]=='S')
		{MAX = (MAX_S-liged); /*cout << endl << "S" <<endl;*/}
		while (!limit_type.eof())
		{
			limit_type.get(atom, sizeof(atom), '\0');	  
			limit_type.get(type, sizeof(type), '\0');
			//cout << type << "*" <<endl;
			if (type[5]=='2' && atom[1]==ch_l[line_l_count][9] && atom[2]==ch_l[line_l_count][10] && atom[3]==ch_l[line_l_count][11] && atom[4]==ch_l[line_l_count][12] && atom[5]==ch_l[line_l_count][13])
			{
				MAX--;
				//cout << "achei uma dupla!!!" << endl; 
			}
			//cout << "type: " << type[4] << type[5] << endl;
		}
		limit_type.clear();              // forget we hit the end of file
		limit_type.seekg(0, ios::beg);   // move to the start of the file
		while (string_line_count<l_TOTAL)
		{
			if ((ch_l[line_l_count][9]!='X')&&(ch_l[line_l_count][9]==ch_l[string_line_count][9])&& (ch_l[line_l_count][10]==ch_l[string_line_count][10])&& (ch_l[line_l_count][11]==ch_l[string_line_count][11])&& (ch_l[line_l_count][12]==ch_l[string_line_count][12])&& (ch_l[line_l_count][13]==ch_l[string_line_count][13]))
			{
				found++; //encontrou um atomo igual ao da string selecionada
				if (found > MAX)
				{
					ch_p[string_line_count][12]='X'; //marca este atomo para excluir da procura
					ch_l[string_line_count][9]='X';
					dist_line[string_line_count][0]='X';
				}
			}
			string_line_count++;
		}
		//cout << endl << "F" << found << "*" << ch_l[line_l_count][9] << endl;
		string_line_count = 0; //reseta as linhas que vao serem comparadas com a selecionada
		found = 0; //reseta nmero de encontrados
		line_l_count++; //seleciona a proxima string a ser comparada com as outras
		//d_list.close();
	}
	line_p_count = 0;
	line_l_count = 0; 
	string_line_count = 0;
	while (string_line_count<l_TOTAL)
	{
		if (ch_l[string_line_count][9]!='X'){refined++;}
		string_line_count++;  
	}
	//polscore.open("polscore.txt", ios::app);
	//cout << "Refined result (HB): " << refined << "\t" << endl; 
	//polscore << refined << endl;
	//polscore.close();
	string_line_count = 0;
}

void deleta_temp()
{
	remove("protein.tmp");
	remove("ligand.tmp");
	remove("bonds.tmp");
	remove("p_result.tmp");
	remove("l_result.tmp");
	remove("P_DIST.tmp");
	remove("L_DIST.tmp");
	remove("P_ROOT_0.tmp");
	remove("L_ROOT_0.tmp");
	remove("P_ROOT.tmp");
	remove("L_ROOT.tmp");
	remove("p_result.tmp");
	remove("l_result.tmp");
	remove("bond_count.txt");
	remove("dist_result.tmp");
	remove("limit_l.tmp");
	remove("limit_p.tmp");
	remove("limit_type.tmp");
	remove("HC_lig_contact.tmp");
	remove("HC_lig_contact_lim.tmp");
	remove("RT_found.tmp");
}
//-----------------------------------------------------------------
//---------------------------------------------------------------
//--------------------------------------------------------------
void calcula_RT() // pronto 
{
	ifstream limit_type;
	limit_type.open("limit_type.tmp");
	char line[8];
	int line_count = 0;
	char type[8];
	int single_bond_count = 0;
	int double_bond_count = 0;
	int ring_count = 0;
	//conta numero de linhas
	while (!limit_type.eof())
	{
		limit_type.get(line, sizeof(line), '\0');
		//cout << "*" << line << "*";
		line_count++;
	}
	char atom_0[line_count][8];
	char atom_1[line_count][8];
	int total_atom_l = line_count;
	//cout << line_count << endl << line_count/2 << endl;
	//reinicia leitura do arquivo
	limit_type.clear();              // forget we hit the end of file
	limit_type.seekg(0, ios::beg);   // move to the start of the file
	int current_line = 0;
	while (!limit_type.eof())
	{
		limit_type.get(atom_1[current_line], sizeof(line), '\0');
		//cout <<  endl << "@" << atom_1[current_line] << "@" ;
		current_line++;
	}
	limit_type.close();
	int lig_line = 0;
	ifstream RT_found;
	RT_found.open("RT_found.tmp");
	char rt_line[8];
	//i = stringObject1.find("Pure");
	// if(i!=string::npos) {
	while (!RT_found.eof())
	{
		RT_found.get(rt_line, sizeof(rt_line), '\0');
		lig_line++; //cout << rt_line << "*" << endl;
	}
	int total_lig_line = lig_line;
	char rt_string[total_lig_line][8];
	lig_line = 0;
	//reinicia leitura do arquivo
	RT_found.clear();              // forget we hit the end of file
	RT_found.seekg(0, ios::beg);   // move to the start of the file
	//pega os numeros dos atomos que estao dentro do permitido para o RT e elimina os repetidos da lista
	while (!RT_found.eof())
	{
		RT_found.get(rt_line, sizeof(rt_line), '\0');
		strcpy(rt_string[lig_line], rt_line); //cout << rt_line << "*";
		lig_line++;
	}
	//deleta repetidos
	int i = 0;
	int j = 0;
	char aux[lig_line][8];
		for(i = 0; i < lig_line; i++)
				strcpy(aux[i],rt_string[i]);  /* copiamos a matriz a aux */
		for(i = 0; i < lig_line; i++)  
		{ for(j = 0; j < lig_line; j++)   
			{if(!strcmp(rt_string[i],aux[j]) && i != j )  /*i != j */
					{strcpy(aux[i],"REP");}
			}    
		}               
	lines_total2 = 0; 
	lines_total2--;     
	for(i = 0, j = 0; i < lig_line; i++)      /*modificam o original*/
	if(strcmp(aux[i],"REP"))    /* se e rep nao copiamos */
	{strcpy(rt_string[j++],aux[i]);lines_total2++;}
	int current_line_0 = 0;
	int current_line_1 = 0;
	lig_line = 0;
	//carrega tipos de atomos
	ifstream ligand;
	ligand.open("ligand.tmp");
	char ligand_line[55];
	char atom_type[line_count/2][12];
	lig_line = 0;
	while (!ligand.eof())
	{
		ligand.get(ligand_line, sizeof(ligand_line), '\0');
		atom_type[lig_line][0] = ligand_line[1]; atom_type[lig_line][1] = ligand_line[2]; atom_type[lig_line][2] = ligand_line[3]; 
		atom_type[lig_line][3] = ligand_line[4]; atom_type[lig_line][4] = ligand_line[5]; atom_type[lig_line][5] = ligand_line[6]; atom_type[lig_line][6] = ligand_line[7]; 
		atom_type[lig_line][7] = ligand_line[47]; atom_type[lig_line][8] = ligand_line[48]; atom_type[lig_line][9] = ligand_line[49];
		atom_type[lig_line][10] = ligand_line[50]; atom_type[lig_line][11] = ligand_line[51];
		lig_line++;
	}
	lig_line = 0;
	ligand.close();
	//ja carregou todos os atomos agora vai comecar a comparar e contar o numero de ligacoes simples.
	int s_bond = 0;
	//abre arquivo mol2 para definir o tipo de cada atomo
	ifstream mol2;
	mol2.open("ligand.tmp");
	int type_l = 0;
	char mol2_line[54];
	while (!mol2.eof())
	{
		mol2.get(mol2_line, sizeof (mol2_line), '\0');
		type_l++;
	}
	//reinicia leitura do arquivo
	mol2.clear();              // forget we hit the end of file
	mol2.seekg(0, ios::beg);   // move to the start of the file
	char mol2_type[type_l][10];
	type_l = 0;
	while (!mol2.eof())
	{
		mol2.get(mol2_line, sizeof (mol2_line), '\0');
		mol2_type[type_l][0] = mol2_line[0];
		mol2_type[type_l][1] = mol2_line[1];
		mol2_type[type_l][2] = mol2_line[2];
		mol2_type[type_l][3] = mol2_line[3];
		mol2_type[type_l][4] = mol2_line[4];
		mol2_type[type_l][5] = mol2_line[5];
		mol2_type[type_l][6] = mol2_line[6];
		mol2_type[type_l][7] = mol2_line[47];
		mol2_type[type_l][8] = mol2_line[48];
		mol2_type[type_l][9] = mol2_line[49];
		type_l++;
	}
	mol2.close();
	//cout << "***************";
	int sp2_count = 0;
	//inserir funcao para encontrar numero do atomo dentro da lista de permitidos RT_string
	//wf: i = rt_string.find(atom_0[current_line_0], 6);
	//if(i==string::npos) 
	//{
	//while(i==string::npos){cout << "not found" << endl; current_line_0++;}
	//}
	//if(i!=string::npos){cout << "found at: " << i << endl;
	int atom_count = 0;
	int marker = 0;
	int pos = 0;
	int pair_count1 = 0;
	int pair_count2 = 1;
	int num_atom1;
	int num_atom2;
	char atom1[7];
	char atom2[7];
	int number_count = 0;
	long pair = 0;
	int next = 0;
	while (current_line_0 < lines_total2)
	{
		//conta quantas ligacoes estao envolvidas com um mesmo atomo
		if (marker < line_count && current_line_0 < lines_total2)
		{
			while (marker < line_count)
			{
				number_count++; //contador de numero de atomos
				next = (marker+1); //posicao onde o tipo de ligacao e descrita
				//se os atomos sao iguais
				if (!strncmp(rt_string[current_line_0], atom_1[marker], 6)){atom_count++;}
				//se os atomos sao iguais e a ligacao eh simples
				if (!strncmp(rt_string[current_line_0], atom_1[marker], 6) && atom_1[next][5]=='1') {
					//se o atomo for igual ao selecionado verifica se os dois atomos sao sp3 ou sp2
					//copia somente o numero dos atomos para uma string
					atom1[0] = atom_1[marker][0]; 
					atom1[1] = atom_1[marker][1];
					atom1[2] = atom_1[marker][2];
					atom1[3] = atom_1[marker][3];
					atom1[4] = atom_1[marker][4];
					atom1[5] = atom_1[marker][5];
					num_atom1 = atoi(atom1);
					if (mol2_type[num_atom1][9]=='2'){sp2_count++;}
					//----------se o segundo atomo eh impar = primeiro de um par
					if (number_count%2 != 0){pair = 2;}
					//----------se o segundo atomo eh par = segundo de um par
					if (number_count%2 == 0){pair = (-2);}
					atom2[0] = atom_1[marker+pair][0]; 
					atom2[1] = atom_1[marker+pair][1];
					atom2[2] = atom_1[marker+pair][2];
					atom2[3] = atom_1[marker+pair][3];
					atom2[4] = atom_1[marker+pair][4];
					atom2[5] = atom_1[marker+pair][5];
					num_atom2 = atoi(atom2);
					if (mol2_type[num_atom2][9]=='2'){sp2_count++;}
					if (sp2_count < 2 && mol2_type[num_atom2][7]!='H' && mol2_type[num_atom1][7]!='H'){s_bond++;} //achou um rotor!!!!
					//	if (sp2_count > 1 && mol2_type[num_atom2][7]!='H') {cout << "2 sp2s";}
				}//se os atomos sao iguais (encontrou um igual ao selecionado na lista)
				marker+=2;
			} 
		}
		marker = 0;
		//-----------------------------------------------------------------
		/*cout << "*" << endl;
		cout << "single bonds" << s_bond << endl;
		cout << "total number of bonds: " << atom_count << endl;
		cout << "number of sp2 atoms: " << sp2_count << endl;
		cout << "*" << endl;
		*/
		if (atom_count > 1)
		{
			if (s_bond == 0) {RT += 0;}
			if (s_bond == 1 || s_bond >= 3){RT += 0.5;}
			if (s_bond == 2) {RT += 1;}
		}
		//cout << "Temporary RT: " << RT <<  " rotors: " << s_bond << " atom_count: " << atom_count << endl;
		s_bond = 0; 
		marker = 0;
		sp2_count = 0;
		atom_count = 0; //zera o numero de atomos econtrados iguais ao selecionado
		current_line_0++;
	}//while (current_line<=line_count)
	//} 
}
void result_score_calc(){ //pronto
	int HB = 0; //HB = hydrogen bonds = C
	HB = refined;
	ofstream result_score;
	char result_file_name[100];
	vector<string> v;
	v.clear();
	split(ligand_name, '.', v); 
	strcpy(result_file_name, v[0].c_str());
	strcat(result_file_name, "_result.txt");
	result_score.open(result_file_name, ios::app);

	//cout << HC_total2 << "\t" << VDW_total << "\t" << RT << "\t" << HB << "\t" << hydrogen_B << "\t" << repulsive << london << endl;
	//cout << "Hydrophobic_contact:         " << HC_total2 << endl;
	//cout << "Van_der_waals:               " << VDW_total << endl;
	//cout << "Deformation_effect:          " << RT << endl;
	//cout << "Hydrogen_bonds_refined_(HB):         " << HB << endl;
	//cout << "Hydrogen_bonds_non_filtered_(HB):         " << hydrogen_B << endl;
	//cout << endl;
	//cout << "---------------------------------------------------------------------------------" << endl;

	result_score << "Hydrophobic contacts:        " << HC_total2 << endl;
	result_score << "Van der waals:               " << VDW_total << endl;
	result_score << "Deformation effect:          " << RT << endl;
	result_score << "Hydrogen bonds (HB):         " << hydrogen_B << endl;
	result_score << "Repulsive VDW score:         " << repulsive << endl;
	result_score << "London dispersion force:     " << london << endl;
	
	double test = 1;
	double test2 = 5;
	double result = summation(test, test2);
	
	//cout << "result: " << endl << result << "char" << protein_name << endl << ligand_pdb << endl;
	double ASA1 = 0;
	double ASA2 = 0;
	string sendme;
	char result1[100];
	char result2[100];
	FILE* fp1;
	FILE* fp2;
	char str1[100];
	strcpy(str1, "python asa/asa.py ");
	char str2[100];
	strcpy(str2, "python asa/asa.py ");
	fp1 = popen(strcat(str1, protein_name), "r");
	fread(result1,1,sizeof(result1),fp1);
	//fclose (fp1);
	pclose(fp1);
	//printf("%s",result2);
	ASA1 = atof(result1);
	fp2 = popen(strcat(str2, ligand_pdb), "r");
	fread(result2,1,sizeof(result2),fp2);
	//fclose (fp2);
	pclose(fp2);
	//printf("%s",result2);
	ASA2 = atof(result2);
	char str3[100];
	strcpy(str3, ligand_name);
	strcat(str3, ".interaction_terms.txt");
	char type[100];
	FILE* fp3;
	char result3[10];
	char str4[100];
	if (strcmp(ligand_type,"DNA") == 0 || strcmp(ligand_type,"RNA") == 0 || strcmp(ligand_type,"DNA/RNA") == 0 || strcmp(ligand_type,"nucleotide") == 0 || strcmp(ligand_type,"protein-DNA") == 0 || strcmp(ligand_type,"protein-RNA") == 0 || strcmp(ligand_type,"protein-DNA/RNA") == 0 || strcmp(ligand_type,"protein-nucleotide") == 0 || strcmp(ligand_type,"nucl") == 0){
		strcpy(str4, "Rscript R/DNA/script.R ");
		strcpy(type,"DNA");
	//cout << "found DNA" << endl;
	}else{
		if (strcmp(ligand_type,"small") == 0 || strcmp(ligand_type,"small molecule") == 0 || strcmp(ligand_type,"drug") == 0 || strcmp(ligand_type,"prtoein-small") == 0 || strcmp(ligand_type,"protein-small molecule") == 0 || strcmp(ligand_type,"protein-drug") == 0){
			strcpy(str4, "Rscript R/small_molecule/script.R ");
			//cout << "found small molecule" << endl;
			strcpy(type,"small");
		}else{
			if (strcmp(ligand_type,"protein") == 0 || strcmp(ligand_type,"protein-protein") == 0 || strcmp(ligand_type,"peptide") == 0 || strcmp(ligand_type,"polypetide") == 0 || strcmp(ligand_type,"prot") == 0){
			strcpy(str4, "Rscript R/protein/script.R ");
			//cout << "found protein" << endl;
			strcpy(type,"protein");
			}else{
				strcpy(str4, "Rscript R/DNA/script.R ");
				if (strcmp(ligand_type,"unknown") == 0){
				cout << "No ligand type provided, using default scoring function (protein-DNA/RNA).\n";
				strcpy(type,"DNA");
				}else{
					cout << "Found unknown ligand type parameter " << ligand_type << ", using default (DNA).\n";
					strcpy(type,"DNA");
				}
			}
		}
	}
	//start: calculates surface tension and hydrophobicity
	string hydro_file = "hydrophobicity.param";
	string tension_file = "tension.param";
	map<string,double> hydro_map, tension_map;
	hydro_map = create_map_tables(hydro_file);
	tension_map = create_map_tables(tension_file);
	//cout << "hydro_ALA " << hydro_map["ALA"] << endl;
	//cout << "tension_ALA " << tension_map["ALA"] << endl;
	string protein_file_name = "protein.tmp";
	string p_dist_file_name = "P_DIST.tmp";
	double total_hydrophobicity = 0;
	double total_surface_tension = 0;
	double contact_hydrophobicity = 0;
	double contact_surface_tension = 0;
	hydrophobicity = 0;
	surface_tension = 0;
	surface_tension_hydrophobicity_calculator(hydro_map, tension_map, protein_file_name);
	total_hydrophobicity = hydrophobicity;
	total_surface_tension = surface_tension;
	//cout << "total hydrophobicity : " << total_hydrophobicity << endl;
	//cout << "total surface tension: " << total_surface_tension << endl;
	hydrophobicity = 0;
	surface_tension = 0;
	surface_tension_hydrophobicity_calculator(hydro_map, tension_map, p_dist_file_name);
	contact_hydrophobicity = hydrophobicity;
	contact_surface_tension = surface_tension;
	//cout << "contact hydrophobicity : " << contact_hydrophobicity << endl;
	//cout << "contact surface tension: " << contact_surface_tension << endl;
	//end: calculates surface tension and hydrophobicity
	ofstream interaction_terms(str3);
	if(strcmp(type,"DNA") == 0){
		//DNA
		interaction_terms << "V2\tV3\tV4\tV5\tV6\tV7\tV18\tV19\tV20\tV21\tV22\tV23"  << endl << HC_total2 << "\t" << VDW_total << "\t" <<  RT << "\t" << hydrogen_B << "\t" << ASA1 << "\t" << ASA2 << "\t" << repulsive << "\t" << london << "\t" << contact_hydrophobicity << "\t" << total_hydrophobicity << "\t" << contact_surface_tension << "\t" << total_surface_tension << endl;
	}
	if(strcmp(type,"small") == 0){
		//small
		interaction_terms << "V2\tV3\tV4\tV5\tV6\tV7\tV17\tV18\tV20\tV21\tV22\tV23"  << endl << HC_total2 << "\t" << VDW_total << "\t" <<  RT << "\t" << hydrogen_B << "\t" << ASA1 << "\t" << ASA2 << "\t" << repulsive << "\t" << london << "\t" << contact_hydrophobicity << "\t" << total_hydrophobicity << "\t" << contact_surface_tension << "\t" << total_surface_tension << endl;
	}
	if(strcmp(type,"protein") == 0){
		//protein
		interaction_terms << "V2\tV3\tV4\tV5\tV6\tV7\tV16\tV17\tV20\tV21\tV22\tV23"  << endl << HC_total2 << "\t" << VDW_total << "\t" <<  RT << "\t" << hydrogen_B << "\t" << ASA1 << "\t" << ASA2 << "\t" << repulsive << "\t" << london << "\t" << contact_hydrophobicity << "\t" << total_hydrophobicity << "\t" << contact_surface_tension << "\t" << total_surface_tension << endl;
	}
	interaction_terms.close();    
	//cout << strcat(str4, str3) << endl;
	fp3 = popen(strcat(str4, str3), "r");
	fread(result3,1,sizeof(result3),fp3);
	pclose(fp3);
	cout << atof(result3) << endl;
	result_score << "predicted pKD:               " << atof(result3) << endl; 	 
}

int main(int argc, char *argv[]) // Pronto
{
    if (argc < 4){
		cout << "usage: GLM-Score <protein_file PDB> <ligand_file PDB> <ligand_file MOL2> <Ligand Type> \n\n<Ligand Type> is the type of ligand molecule, which supports 3 alternative options:\n\nDNA (or RNA), small-molecule (or drug), and protein.\nIf no <Ligand Type> is provided, protein-DNA/RNA scoring function will be used as default.\n\n";
		exit(1);
    }
	strcpy(protein_name, argv[1]);
	strcpy(ligand_name, argv[3]);
	strcpy(ligand_pdb, argv[2]);
	if (argc == 4){
		strcpy(ligand_type, "unknown");
	}
	if (argc == 5){
		strcpy(ligand_type, argv[4]);
	}
	// 1) processa entreadas
	// [JC] Open the file ligand.mol2
	open_ligand();
	// [JC] Open PDF file
	open_protein();

	// 2) Calculos
	// 2.1 [JC] Hydrophobic contacts (HCs) / van der Waals interactions (VDWs) / Repulsive interactions (RIs)	
	distancia();

	// 2.2 [JC]
	raiz_proteina();

	raiz_ligante();
	P_centro_geom();
	L_centro_geom();
	angulos();
	// 3) Saida
	salva_ligante();
	saida_PDB();
	salva_proteina();
	// 4) refinamento do HB
	if (hydrogen_B!=0)
	{
		limit();
	}
	calcula_RT();
	result_score_calc();
	deleta_temp();
	//system ("pause");
}
