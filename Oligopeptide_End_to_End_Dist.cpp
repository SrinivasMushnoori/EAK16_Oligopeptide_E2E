#include<iostream>//cout
#include<string>//string
#include<vector>//vector
#define _USE_MATH_DEFINES
#include<math.h>//fabs
#include<fstream>//ofstream
#include<algorithm>//find
using namespace std;

void readump01(string file_name, vector<int> &steps, int &N, vector<double> &xlo, vector<double> &xhi, vector<double> &ylo, vector<double> &yhi, vector<double> &zlo, vector<double> &zhi, vector<vector<int> > &id, vector<vector<int> > &mol, vector<vector<int> > &type, vector<vector<double> > &q, vector<vector<double> > &x, vector<vector<double> > &y, vector<vector<double> > &z, vector<vector<double> > &ix, vector<vector<double> > &iy, vector<vector<double> > &iz, vector<vector<double> > &masses, vector<vector<double> > &vx, vector<vector<double> > &vy, vector<vector<double> > &vz, vector<vector<double> > &fx, vector<vector<double> > &fy, vector<vector<double> > &fz);
//double periodic(double init,double fin,double lo,double hi);
void readata01(string data_file_name, vector<int> &bondIDs, vector<int> &bond_types, vector<vector<int> > &bonds, vector<int> &angleIDs, vector<int> &angle_types, vector<vector<int> > &angles, vector<int> &dihedralIDs, vector<int> &dihedral_types, vector<vector<int> > &dihedrals);

int main()
{
	//VARIABLE DEFINITIONS

	int N, files;
	vector<string> file_names;
	vector<int> steps, tmp_steps, bond_types, bondIDs, angle_types, angleIDs, dihedralIDs, dihedral_types;
	vector<double> xlo, xhi, ylo, yhi, zlo, zhi, tmp_xlo, tmp_xhi, tmp_ylo, tmp_yhi, tmp_zlo, tmp_zhi;
	vector<vector<int> > id, mol, type, tmp_id, tmp_mol, tmp_type;
	vector<vector<int> > bonds, angles, dihedrals;
	vector<vector<double> > q, x, y, z, ix, iy, iz, masses, vx, vy, vz, fx, fy, fz, tmp_q, tmp_x, tmp_y, tmp_z, tmp_ix, tmp_iy, tmp_iz, tmp_masses, tmp_vx, tmp_vy, tmp_vz, tmp_fx, tmp_fy, tmp_fz;
	vector<double> mediancalc;
	int  backbone_bond_types, back_length;
	string data_file_name;
	double dx, dy, dz;
	vector<int> backbone_bonds;
	vector<vector<double> > E2ED, E2E;
	vector<double> E2E_avg;
	//END OF VARIABLE DEFINITIONS


	//USER INPUT DUMP/DATA FILE 
	cout << "number of dump files = ";
	cin >> files;
	file_names.resize(files);
	for (int file = 0; file < files; file++)
	{
		cout << "name of dump file " << file + 1 << ": ";
		cin >> file_names[file];
	}


	cout << "Data file name: ";
	cin >> data_file_name;

	cout << "peptide backbone length = ";
	cin >> back_length;

	//cout << "number of bond types in peptide backbone = ";
	//cin >> backbone_bond_types;
	//for (int bond_counter = 0; bond_counter < backbone_bond_types; bond_counter++)
	//{
	//	backbone_bonds.push_back(0); //Vector called backbone_bonds contains the types of relevant bonds
	//	cout << "Backbone bond type " << bond_counter + 1 << " = ";
	//	cin >> backbone_bonds[bond_counter];
	//}



	//ASSEMBLING ONE MASSIVE DUMPFILE

	readump01(file_names[0], steps, N, xlo, xhi, ylo, yhi, zlo, zhi, id, mol, type, q, x, y, z, ix, iy, iz, masses, vx, vy, vz, fx, fy, fz);

	for (int file = 1; file < files; file++)
	{
		readump01(file_names[file], tmp_steps, N, tmp_xlo, tmp_xhi, tmp_ylo, tmp_yhi, tmp_zlo, tmp_zhi, tmp_id, tmp_mol, tmp_type, tmp_q, tmp_x, tmp_y, tmp_z, tmp_ix, tmp_iy, tmp_iz, tmp_masses, tmp_vx, tmp_vy, tmp_vz, tmp_fx, tmp_fy, tmp_fz);
		steps.insert(steps.end(), tmp_steps.begin() + 1, tmp_steps.end());
		xlo.insert(xlo.end(), tmp_xlo.begin() + 1, tmp_xlo.end());
		xhi.insert(xhi.end(), tmp_xhi.begin() + 1, tmp_xhi.end());
		ylo.insert(ylo.end(), tmp_ylo.begin() + 1, tmp_ylo.end());
		yhi.insert(yhi.end(), tmp_yhi.begin() + 1, tmp_yhi.end());
		zlo.insert(zlo.end(), tmp_zlo.begin() + 1, tmp_zlo.end());
		zhi.insert(zhi.end(), tmp_zhi.begin() + 1, tmp_zhi.end());
		id.insert(id.end(), tmp_id.begin() + 1, tmp_id.end());
		mol.insert(mol.end(), tmp_mol.begin() + 1, tmp_mol.end());
		type.insert(type.end(), tmp_type.begin() + 1, tmp_type.end());
		q.insert(q.end(), tmp_q.begin() + 1, tmp_q.end());
		x.insert(x.end(), tmp_x.begin() + 1, tmp_x.end());
		y.insert(y.end(), tmp_y.begin() + 1, tmp_y.end());
		z.insert(z.end(), tmp_z.begin() + 1, tmp_z.end());
		ix.insert(ix.end(), tmp_ix.begin() + 1, tmp_ix.end());
		iy.insert(iy.end(), tmp_iy.begin() + 1, tmp_iy.end());
		iz.insert(iz.end(), tmp_iz.begin() + 1, tmp_iz.end());
		masses.insert(masses.end(), tmp_masses.begin() + 1, tmp_masses.end());
		vx.insert(vx.end(), tmp_vx.begin() + 1, tmp_vx.end());
		vy.insert(vy.end(), tmp_vy.begin() + 1, tmp_vy.end());
		vz.insert(vz.end(), tmp_vz.begin() + 1, tmp_vz.end());
		fx.insert(fx.end(), tmp_fx.begin() + 1, tmp_fx.end());
		fy.insert(fy.end(), tmp_fy.begin() + 1, tmp_fy.end());
		fz.insert(fz.end(), tmp_fz.begin() + 1, tmp_fz.end());
	}

	//READING THE DATA FILE, SCREENING
	readata01(data_file_name, bondIDs, bond_types, bonds, angleIDs, angle_types, angles, dihedralIDs, dihedral_types, dihedrals);
	vector<int> filtered_bonds1, filtered_bonds2, pep_start, pep_end;
	vector< vector<int> > peptides;

	//cout<<"\nSize of bond types :: "<<bond_types.size()<<" versus "<< bondIDs.size()<<"\n";

	//filtered_bonds1.push_back(0);
	//filtered_bonds2.push_back(0);


	for (int bond_count = 0; bond_count < bondIDs.size(); bond_count++)
	{
		
		if (bond_types[bond_count] == 1)
		{
			//cout << "Bond count condition satisfied , bond_count= "<< bond_count<<", bondIDs[bond_count]= "<<bondIDs[bond_count] << ", bond_types[bond_count]= "<< bond_types[bond_count]<<"\n";
			filtered_bonds1.push_back(bonds[bond_count][0]);
			filtered_bonds2.push_back(bonds[bond_count][1]);
			//filtered_bonds1 and 2 are basically column 1 and 2 of the "bonds" vector
		}
	}



		
	//Using the fact that the number of beads/peptide is known, and the fact that peptides are 
	//sequential, we know that if A is the atom ID of the starting peptide bead, A+15 will be the corresponding terminal bead

	sort(filtered_bonds1.begin(), filtered_bonds1.end());
	sort(filtered_bonds2.begin(), filtered_bonds2.end());



	for (int i = 0; i < filtered_bonds1.size(); i++)
	{
		//cout << filtered_bonds1[i] << " " << filtered_bonds2[i] << endl;
	}



	set_difference(filtered_bonds1.begin(), filtered_bonds1.end(), filtered_bonds2.begin(), filtered_bonds2.end(), inserter(pep_start, pep_start.begin()));
	int len = pep_start.size();
	cout << "pep_start length is " << len << endl;


	for (int i = 0; i < len; i++)
	{

		//cout << "loop triggered" << endl;
		cout << "pep_start[k] is " << pep_start[i] << endl;

	}




	


	E2E.resize(steps.size());
	E2E_avg.resize(steps.size());
	vector<double> Std_Dev;
	Std_Dev.resize(steps.size());
	//double avdist;
	//double sumdist = 0;
	//int time = 50;
	/*---------------------------------------------------------------------------------------------------------------------
	for (int time = 0; time < steps.size(); time++)
	{

		for (int k = 0; k < pep_start.size(); k++)
		{
			//cout << " check 1: getting this far " << endl;
			E2E[time].resize(pep_start.size());
			int w = pep_start[k] - 1;
			//cout << " check 2: reading pepstart k correctly " << endl;
			int u = w + back_length - 1;


			double xpb = x[time][w] + (ix[time][w] * (xhi[time] - xlo[time])); //beginning bead coordinates
			double ypb = y[time][w] + (iy[time][w] * (yhi[time] - ylo[time]));
			double zpb = z[time][w] + (iz[time][w] * (zhi[time] - zlo[time]));
			//let's see if it reads w for 16
			cout << " w=16 triggercheck " << endl;
			double xpe = x[time][u] + (ix[time][u] * (xhi[time] - xlo[time])); //ending bead coordinates
			double ype = y[time][u] + (iy[time][u] * (yhi[time] - ylo[time]));
			double zpe = z[time][u] + (iz[time][u] * (zhi[time] - zlo[time]));

			double xp, yp, zp;
			xp = xpe - xpb;
			yp = ype - ypb;
			zp = zpb - zpe;
			double End_to_End = sqrt((xp*xp) + (yp*yp) + (zp*zp));

			//cout << " check 3: calculations happening " << endl;
			E2E[time][k] = End_to_End;
			//E2E[time].push_back(End_to_End);
			sumdist = sumdist + E2E[time][k];
			cout << " check 3: sumdist being written correctly " << endl;
			vector<double> avdist[time] = sumdist / pep_start.size();

			//cout << " End to end dist for chain " << k << " is " << End_to_End << endl;
			//cout << " Current average is " << avdist << endl;
			//cout << " loop count is  " << k << endl;


		}

	}

	//cout << " Average E2Ed for current timestep is " << avdist << endl;

        double tempvar=E2E[0][0];  
        cout << "E2E dist for first chain (starting at 84022) at the first timestep is " << tempvar;



	-------------------------------------------------------------------------------------------------------------------*/



	
	//DISTANCE CALCULATED IS TOO SMALL
	for (int t = 0; t < steps.size(); t++)
	//for (int k = 0; k < pep_start.size(); k++)
	{
		double sum = 0;
		double squarestemp = 0;
		//for (int t = 0; t < steps.size(); t++)
		for (int k = 0; k < pep_start.size(); k++)
		{
			
			E2E[t].resize(pep_start.size());
			int l = pep_start[k]-1;
			int m = l + back_length - 1;
			//cout << l << endl;



			dx = x[t][l] + (ix[t][l] * (xhi[t] - xlo[t])) - x[t][m] - (ix[t][m] * (xhi[t] - xlo[t]));
			dy = y[t][l] + (iy[t][l] * (yhi[t] - ylo[t])) - y[t][m] - (iy[t][m] * (yhi[t] - ylo[t]));
			dz = z[t][l] + (iz[t][l] * (zhi[t] - zlo[t])) - z[t][m] - (iz[t][m] * (zhi[t] - zlo[t]));

			double e2e_d = sqrt((dx*dx) + (dy*dy) + (dz*dz));
			E2E[t][k] = e2e_d;
			mediancalc.push_back(E2E[t][k]);
			//cout << y[t][l] << " " << y[t][l+15] << " " << dy << endl; 
			//cout << l << " " << l+15 << endl;   

			//cout << dx << " " << dy << " " << dz << " " << e2e_d << endl;

			sum += e2e_d;
			//cout << e2e_d << " " << sum << endl;
			//cout << e2e_d << endl;
			//squarestemp = squarestemp + (e2e_d - (sum / pep_start.size()))*(e2e_d - (sum / pep_start.size()));
			
			

			E2E_avg[t] = sum / pep_start.size();
			
                        




		}

		//Std_Dev[t] = squarestemp / steps.size();
	}


	sort(mediancalc.begin(), mediancalc.end());

        //for (int i=0;i<mediancalc.size();i++){cout << mediancalc[i]<<endl;}
             


	int mediancalclength=mediancalc.size();
         
        //cout << "mediancalclength/2 remainder is" << mediancalclength%2 << endl;    
        //cout << mediancalc[mediancalclength/2] << " " << mediancalc[1+mediancalclength/2] << endl;
        double midpoints1=mediancalc[mediancalclength/2];
        double midpoints2=mediancalc[1+mediancalclength/2];
        
	if (mediancalclength%2 == 1)
		{
		double median1=mediancalc[(((mediancalclength-1)/2)+1)];

		cout<<"median is "<<median1<<endl;
		}
	else
		{
		double median1=(midpoints1+midpoints2)/2;
		cout<<"median is "<<median1<<endl;
		}



	



        double tempvar=E2E[0][0];  
        cout << "E2E dist for first chain (starting at 84022) at the first timestep is " << tempvar << endl;
	



	double diffsquared = 0;
	
	
	for (int u = 0; u < steps.size(); u++)
	{
		
		for (int v = 0; v < pep_start.size(); v++)

		{
			double mean = E2E_avg[u];
			//cout << E2E_avg[u] << endl;
			double dist = E2E[u][v];
                        //cout << " mean - dist = " << mean - dist << endl;
			double diffsquared; 
                        diffsquared = diffsquared + (mean - dist)*(mean - dist);
                        //cout << diffsquared << endl;
  
               		Std_Dev[u] = sqrt(diffsquared / pep_start.size());
                        diffsquared=0;
			
		}
               
	}
//------------------------------------HISTOGRAM PLOTTER SECTION------------------------------------------//


vector<int> hist;

hist.resize(16);

hist[0]=0;
hist[1]=0;
hist[2]=0;
hist[3]=0;
hist[4]=0;
hist[5]=0;
hist[6]=0;
hist[7]=0;
hist[8]=0;
hist[9]=0;
hist[10]=0;
hist[11]=0;
hist[12]=0;
hist[13]=0;
hist[14]=0;
hist[15]=0;



for (int u = 0; u < steps.size(); u++)
	{
		
		for (int v = 0; v < pep_start.size(); v++)

		{
			
			double dist = E2E[u][v];
                        if (dist<0.499999)
                        {
                        	hist[0]=hist[0]+1;
                        }
                        else if ((dist>0.5) && (dist<0.99999))
                        {
                        	hist[1]=hist[1]+1;
                        }
                        else if ((dist>1) && (dist<1.499999))
                        {
                        	hist[2]=hist[2]+1;
                        }
                        else if ((dist>1.5) && (dist<1.99999))
                        {
                        	hist[3]=hist[3]+1;
                        }
                        else if ((dist>2) && (dist<2.499999))
                        {
                        	hist[4]=hist[4]+1;
                        }
                        else if ((dist>2.5) && (dist<2.99999))
                        {
                        	hist[5]=hist[5]+1;
                        }
                        else if ((dist>3) && (dist<3.499999))
                        {
                        	hist[6]=hist[6]+1;
                        }
                        else if ((dist>3.5) && (dist<3.99999))
                        {
                        	hist[7]=hist[7]+1;
                        }
                        else if ((dist>4) && (dist<4.499999))
                        {
                        	hist[8]=hist[8]+1;
                        }
                        else if ((dist>4.5) && (dist<5.99999))
                        {
                        	hist[9]=hist[9]+1;
                        }
                        else if ((dist>6) && (dist<6.499999))
                        {
                        	hist[10]=hist[10]+1;
                        }
                        else if ((dist>6.5) && (dist<6.99999))
                        {
                        	hist[11]=hist[11]+1;
                        }
                        else if ((dist>7) && (dist<7.499999))
                        {
                        	hist[12]=hist[12]+1;
                        }
                        else if ((dist>7.5) && (dist<7.99999))
                        {
                        	hist[13]=hist[13]+1;
                        }
                        else if ((dist>8) && (dist<8.499999))
                        {
                        	hist[14]=hist[14]+1;
                        }   
                        else if ((dist>8.5) && (dist<8.99999))
                        {
                        	hist[15]=hist[15]+1;
                        }   
                        else if ((dist>9) && (dist<9.499999))
                        {
                        	hist[16]=hist[16]+1;
                        }   
                        else if ((dist>9.5) && (dist<9.99999))
                        {
                        	hist[17]=hist[17]+1;
                        }   
                        else if ((dist>10) && (dist<10.499999))
                        {
                        	hist[18]=hist[18]+1;
                        }   
                        else if ((dist>10.5) && (dist<10.99999))
                        {
                        	hist[19]=hist[19]+1;
                        }   
                        else if (dist>11)
                        {
                        	hist[20]=hist[20]+1;
                        }
                        		
		}
               
	}

for (int i=0;i<16;i++)
{
	cout << hist[i] << endl;
}



//-----------------------------------END OF HISTOGRAM PLOTTER----------------------------------------------//	


	//cout << Std_Dev[45] << endl;

	ofstream E2E_file;
	E2E_file.open("EndToEndDistance.txt", ios::trunc);
	E2E_file << "Steps" << "\t" << "EtoEDist" << "\t" << "stdev" << "\n";

	for (int a = 0; a<steps.size(); a++){
		E2E_file << a << "\t" << E2E_avg[a] << "\t" << Std_Dev[a] << "\n";
	}

	E2E_file.close();
	cout << "press any key & enter to close";
	cin >> file_names[0];
	return 0;
}




















