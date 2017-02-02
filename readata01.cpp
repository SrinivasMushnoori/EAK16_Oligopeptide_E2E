#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include<iostream>//string
#include<vector>//vector
#include<fstream>//ifstream
#include<sstream>//istringstream
using namespace std;

void readata01(string data_file_name,vector<int> &bondIDs,vector<int> &bond_types,vector<vector<int> > &bonds,vector<int> &angleIDs,vector<int> &angle_types,vector<vector<int> > &angles,vector<int> &dihedralIDs,vector<int> &dihedral_types,vector<vector<int> > &dihedrals){
	bondIDs.clear();
	bond_types.clear();
	bonds.clear();
	angleIDs.clear();
	angle_types.clear();
	angles.clear();
	dihedralIDs.clear();
	dihedral_types.clear();
	dihedrals.clear();
	ifstream data;
	data.open(data_file_name.c_str());
	string line;
	while(getline(data,line)){
		istringstream iss(line);
    	string sub;
    	iss>>sub;
		if(sub=="Bonds"){
			getline(data,line);
			getline(data,line);
			while(!line.empty()){
				bondIDs.push_back(0);
				bond_types.push_back(0);
				bonds.push_back(vector<int>(2,0));
				sscanf(line.c_str(),"%i %i %i %i",&bondIDs.back(),&bond_types.back(),&bonds.back()[0],&bonds.back()[1]);
				getline(data,line);
				istringstream iss(line);
				string sub;
				iss>>sub;
			}
		}if(sub=="Angles"){
			getline(data,line);
			getline(data,line);
			while(!line.empty()){
				angleIDs.push_back(0);
				angle_types.push_back(0);
				angles.push_back(vector<int>(3,0));
				sscanf(line.c_str(),"%i %i %i %i %i",&angleIDs.back(),&angle_types.back(),&angles.back()[0],&angles.back()[1],&angles.back()[2]);
				getline(data,line);
				istringstream iss(line);
				string sub;
				iss>>sub;
			}
		}if(sub=="Dihedrals"){
			getline(data,line);
			getline(data,line);
			while(!line.empty()){
				dihedralIDs.push_back(0);
				dihedral_types.push_back(0);
				dihedrals.push_back(vector<int>(4,0));
				sscanf(line.c_str(),"%i %i %i %i %i %i",&dihedralIDs.back(),&dihedral_types.back(),&dihedrals.back()[0],&dihedrals.back()[1],&dihedrals.back()[2],&dihedrals.back()[3]);
				getline(data,line);
				istringstream iss(line);
				string sub;
				iss>>sub;
			}
		}
	}data.close();
}
/*
int main(){
	string file_name;
	vector<int> bondIDs,bond_types,angleIDs,angle_types,dihedralIDs,dihedral_types;
	vector<vector<int> > bonds,angles,dihedrals;
	cout<<"data file: ";
	cin>>file_name;
	readata01(file_name,bondIDs,bond_types,bonds,angleIDs,angle_types,angles,dihedralIDs,dihedral_types,dihedrals);
	for(int i=0;i<bondIDs.size();i++){
		cout<<bondIDs[i]<<"\t"<<bond_types[i]<<"\t"<<bonds[i][0]<<"\t"<<bonds[i][1]<<"\n";
	}for(int i=0;i<angleIDs.size();i++){
		cout<<angleIDs[i]<<"\t"<<angle_types[i]<<"\t"<<angles[i][0]<<"\t"<<angles[i][1]<<"\t"<<angles[i][2]<<"\n";
	}for(int i=0;i<dihedralIDs.size();i++){
		cout<<dihedralIDs[i]<<"\t"<<dihedral_types[i]<<"\t"<<dihedrals[i][0]<<"\t"<<dihedrals[i][1]<<"\t"<<dihedrals[i][2]<<"\t"<<dihedrals[i][3]<<"\n";
	}cout<<"press any key & enter to close";
	cin>>file_name;
	return 0;
}*/
