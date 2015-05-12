// #!/usr/bin/env python
// """
// comes from Trj_modify is used to fit the solute group in trajectory to the center.
// version 1.0 at 14.11.03
// @change:
//     - 2011.12.26.
//         - Add those words.
//         - Modified it.
//     - 2014.06.23.
//         - Update the notes.
// """

#include "string_operate.h"
#include "read_ndx.h"
#include <math.h>
#include <algorithm>

extern "C"
{
#include "xdrfile_xtc.h"
}


void Print_usage(){
    cout<<"Usage: traj-warp trjfile index.ndx output.xtc"<< endl;
}

void Move_2_center(char * trj_file, char * index_file, char * trjout_file)
{

    // map<int,atom> atom_list    =read_pdb_to_atom(top_file);
    int centering_group;
    // atom_list        =copy.deepcopy(Atoms)
    vector <Index_class> index_list       =Read_index_to_Inclass(index_file);
    Print_Index(index_list);
    cout << "Choosing the group for centering:" ;
    cin >> centering_group; 
    vector<int> centering_list  =index_list[centering_group].group_list ;
    int centering_atoms = centering_list.size();

    int solvent_group;
    cout << "Choosing the solvent group:" ;
    cin >> solvent_group; 
    vector<int>solvent_list = index_list[solvent_group].group_list;
    int solvent_atoms = solvent_list.size();

    int natoms,step;
    // int OUTPUT_ATOMS;
    // float OFF;
    float time_temp;
    float p;
    // vector<double> result;
    matrix box;
    rvec *X_in;
    rvec *X_out;
    XDRFILE *xtc_in;
    XDRFILE *xtc_out;
    xtc_in=xdrfile_open(trj_file,"r");
    int read_return=read_xtc_natoms(trj_file,&natoms);
    X_in=(rvec * )calloc(natoms,sizeof(X_in[0]));
    xtc_out = xdrfile_open(trjout_file,"w");
    X_out=(rvec * )calloc(natoms,sizeof(X_in[0]));

    while(1)
    {
        read_return=read_xtc(xtc_in,natoms,&step,&time_temp,box,X_in,&p);
        if(read_return!=0)
        {
            break;
        }
        //if(step % 100000==0)
        //{
        	cout << "step: "<< step << endl;
        //}
        for (int i = 0; i < 3; ++i)
        {
        	for (int j = 0; j < 3; ++j)
        	{
        		cout<<"\t"<<box[i][j]<<"\t";
        	}
        	cout<<endl;
        }

        // calculate the distance between the center of centering_group and center of box [ref_com];
        float ref_com[3];
        for(int i=0;i< centering_atoms;i++)
        {
            for(int j=0;j<3;j++)
            {
                ref_com[j] += X_in[centering_list[i]-1][j];
            }
        }

        for(int i=0;i< 3;i++)
        {
            ref_com[i] = ref_com[i]/centering_atoms - box[i][i]/2;
        }

        // cout <<   
        // ref_com = ref_com - np.array([dimensions[0]/2, dimensions[1]/2, dimensions[2]/2])  
        /*for (int i=0;i <natoms ; i++)
        {

            for(int j=0;j<3;j++)
            {
                X_out[i][j] = X_in[i][j] - ref_com[j];  

                while( X_out[i][j] <0)
                {
                	X_out[i][j] += box[j][j];
                }

                if (X_out[i][j] > box[j][j])
                {
                    X_out[i][j]=fmod(X_out[i][j],box[j][j]);
                }

            }
        }
		*/



        for (int i=0;i <natoms ; i++)
        {

            for(int j=0;j<3;j++)
            {
                X_out[i][j] = X_in[i][j] - ref_com[j];
            }

            //check whether the selected atom is in the centering or solvent group or not !!!
            int atomSerial = i + 1;
            
            
            bool pionterCl=false;
            bool pionterSl=false;
            /*
            // cout<<"pionterSl=\t"<<pionterSl<<"\t pionterCl=\t"<<pionterCl<<endl;
            vector<int>::iterator icl;
			icl = find(centering_list.begin(), centering_list.end(), atomSerial);
			if(icl != centering_list.end())
			{
				pionterCl=true;
			}

			vector<int>::iterator isl;
			icl = find(solvent_list.begin(), solvent_list.end(), atomSerial);
			if(isl != solvent_list.end())
			{
				pionterSl=true;
			}
			
			if (pionterCl or pionterSl)
			{
				continue;
			}

			cout<<"hello!!"<<endl;

			*/

            for (int m = 0; m < centering_atoms; ++m)
            {
                if (atomSerial==centering_list[m])
                    {
                        pionterCl=true;
                        break;
                    }    
            }

            for (int m = 0; m < solvent_atoms; ++m)
            {
                if (atomSerial==solvent_list[m])
                    {
                        pionterSl=true;
                        break;
                    }    
            }

            if (pionterCl or pionterSl)
            {
                continue;
            }

            // cout<<".... i= \t"<<i<<endl;
            for(int j=0;j<3;j++)
            {
                
                while( X_out[i][j] < 0)
                {
                	X_out[i][j] += box[j][j];
                }


                
                if (X_out[i][j] > box[j][j])
                {
                    X_out[i][j]=fmod(X_out[i][j],box[j][j]);
                }


            }

              
        }

        //apply PBC to solvent for each frame;
        cout<<"apply PBC to solvent for each frame"<<endl;
       for (int i = 0; i < solvent_atoms; ++i)
       {
       		if(i % 3 == 0)
       		{
       			int atomIndex = solvent_list[i] - 1;
       			for(int j=0;j<3;j++)
            	{
                	while(X_out[atomIndex][j] < 0)
                	{
                		X_out[atomIndex][j] += box[j][j];
                		X_out[atomIndex+1][j] += box[j][j];
                		X_out[atomIndex+2][j] += box[j][j];

               		}
               		while(X_out[atomIndex][j] > box[j][j])
                	{
                		X_out[atomIndex][j] -= box[j][j];
                		X_out[atomIndex+1][j] -= box[j][j];
                		X_out[atomIndex+2][j] -= box[j][j];

               		}

            	}

       		}
       }
        read_return = write_xtc(xtc_out, natoms, step, time_temp, box, X_out, p);
    }
}


// if __name__=="__main__":
//     if len(sys.argv)!=5:
//         Usage()
//         sys.exit()

//     topol =sys.argv[1] #PRMpbc
//     intrj =sys.argv[2] #TRJpbc_bz2
//     index = sys.argv[3]
//     outtrj = sys.argv[4]
//     Move_2_center(topol,intrj,index,outtrj)    



int main(int argc,char * argv[])
{
    // unit is nm
    char * coor_file;
    char * traj_file;
    char * index_file;
    // char * data_file;
    char * trjout_file;
    float cutoff = 1.0;

    switch(argc)
    {
        case 4:
            // coor_file = argv[1];
            traj_file = argv[1];
            index_file = argv[2];
            trjout_file = argv[3];
            // cutoff = atoi(argv[5]);

            Move_2_center(traj_file,index_file,trjout_file);
            break;

        case 2:
            Print_usage();
            exit(0);

        default:
            Print_usage();
            exit(0);
    }
}

