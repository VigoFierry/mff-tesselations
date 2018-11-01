#define _CRT_SECURE_NO_DEPRECATE

#pragma warning(disable : 4996)


// #include "C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/src/voro++.hh"

#include "Header.h"


using namespace voro; // kvuli zavedeni trid (napr kvuli pouziti slova pre_container)



int main()
{


            int g,h, i;

			long n_add = 0; // citace poctu pridanych/smazanych/posunutych castic (celkovy pocet castic lze odtud pak dopocitat)
			long n_del = 0;
			long n_mov = 0;
                        //long n_chr = 0;
			long n_part;
			long n;

			//int ijk, q;
			std::vector<int> fid;
            int nx, ny, nz;

			pre_container_poly pconp(0, 1, 0, 1, 0, 1, true, true, true);  // true = periodic in given coordinate

			// Krill data window:
			//pre_container_poly pconp(160, 646, 140, 669, 116, 801, true, true, true);  // true = periodic in given coordinate
			//pre_container_poly pconp(160, 646, 140, 669, 116, 801, false, false, false);

			//pconp.import("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/SimulationsPRPAP/study7/src1.0/datacon.txt");
			//pconp.import("C://Users/Vigo/Dokumenty/Visual Studio 2015/Projects/Gibbs-Voronoi/Data/data_L8000div.txt");
			//            pconp.import("../data/Lag2.txt");
            //            pconp.import("../data/dataL_random.txt");
						pconp.import("../Data/dataconHIST.txt");

			//pconp.import("C:/Random tesselations/Data/data Laguerre/sample6/datacon_c.txt");
			//pconp.import("C:/Users/Vigo/Documents/Visual Studio 2015/Projects/Gibbs-Voronoi/Gibbs-Voronoi/dataconHIST.txt");
			// Krill data:
			//pconp.import("C:/Random tesselations/Data/data Laguerre/Krill data/datacon_c.txt");
			pconp.guess_optimal(nx, ny, nz);  // guess
			// Set up the container class and import the particles from the pre-container

			// original container >
			container_poly conp(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8);
			// copy of the container >
			container_poly conpc(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8);

			// Krill data:
			//container_poly conp(160, 646, 140, 669, 116, 801, nx, ny, nz, true, true, true, 8);
			//container_poly conp(160, 646, 140, 669, 116, 801, nx, ny, nz, false, false, false, 8);

			pconp.setup(conp);		// import do originalu
			pconp.setup(conpc);		// import do kopie
			// timto jsou vytvoreny dve stejne kopie containeru, nadale je potreba udrzovat tyto kopie identicke

			n = conp.total_particles();  // indikator maximalniho pouziteho ID (vcetne vektoru fid)


			// PARAMETERS
			std::vector<double> hpar;
			hpar.resize(3);
			//hpar.resize(4);
			hpar[0] = 0.002;	// alfa
			hpar[1] = 0.08;	// beta
			hpar[2] = 0;	// B
			//hpar[3] = 0;	// iota		(optional)
			double alfa_e, beta_e, B_e, iota_e, theta_e, zet_e; // promenne pro odhady parametru

			// unconstrained program:	alfa = 0; alfa_e = alfa;	beta = 3000; beta_e = beta;		B = 0; B_e = 0;

			//alfa = 0.02;  // "alfa = h_min"
			//beta = 0.05;	  // "beta = h_max"
			//B = 1.5;    // "B = h_max^3 / vol"

			std::vector<double> theta;
			theta.resize(3);
			theta[0] = 1;
			theta[1] = 1;
			theta[2] = 1;
			// pro theta zaporne musime davat pozor na STABILITU
			//double zet = 0.000001;
			double zet = 1000000000;
			double sigma = pow(0.015, 2);

			n_part = n;

			//std::cout << "Empty cells: " << empty_cells(conp) << " " << empty_cells(conpc) << "\n";
			//std::cout << "No of boxes: " << conp.nxyz << " " << conpc.nxyz << "\n";
			//delete_empty(conp);
			//std::cout << " next \n";
			//delete_empty(conpc);
			//std::cout << "Empty cells: " << empty_cells(conp) << " " << empty_cells(conpc) << "\n";
			//std::cout << "No of boxes: " << conp.nxyz << " " << conpc.nxyz << "\n";

			
			std::cout << "Total: " << conp.total_particles() << " \n";
			std::cout << "Empty: " << empty_cells(conp) << " \n";
			std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

			delete_empty(conp); delete_empty(conpc);

			std::cout << "Total: " << conp.total_particles() << " \n";
			std::cout << "Empty: " << empty_cells(conp) << " \n";
			std::cout << "Non-empty: " << nonempty_cells(conp) << " \n";

			if (feasibility(conp, hpar[0], hpar[1], hpar[2])) { std::cout << " Feasible \n"; }
			else { std::cout << " Unfeasible \n"; }

			hardcore_estim(conp, alfa_e, beta_e, B_e);
			std::cout << alfa_e << " " << beta_e << " " << B_e << " \n";

			if (feasibility(conp, alfa_e, beta_e, B_e)) { std::cout << " Feasible under estimates \n"; }
			else { std::cout << " Unfeasible under estimates \n"; }

			std::vector<double> hpar_estim;
			hpar_estim.resize(3);
			hpar_estim[0] = alfa_e;	// alfa
			//hpar_estim[0] = 0;	// alfa unspecified
			hpar_estim[1] = beta_e;	// beta
			hpar_estim[2] = B_e;	// B
			//hpar_estim[2] = 0;	// B unspecified

			std::vector<double> vec;  // removable coefficients


			

		// vypocet vstupni informace (je-li potreba, jinak vynechat) - vyberove prumery, histogramy charakteristik
			con_info info;
			info.get_tp(conp);
			//std::cout << info.tp_bef << "\n";
			info.get_meansum(conp);
			//std::cout << info.mean_bef.size() << " " << info.mean_aft.size() << "\n";
			//std::cout << info.mean_bef[0] << " " << info.mean_aft[0] << "\n";
			//info.get_varsum(conp);
			//std::cout << "var " << info.var_bef[0] / (info.tp_bef - 1) << " \n";

			histogram hist_vol, hist_nof;
			hist_nof.create_hist_int(conp); // nof
			hist_vol.create_hist_double(conp); // volume
			//std::cout << hist.sp << " " << hist.step << " " << hist.oc.size() << " " << hist.oc[0] << " " << hist.ocm << "\n";
			info.hist2_bef = hist_vol;
			info.hist2_aft = hist_vol;
			info.hist_bef = hist_nof;
			info.hist_aft = hist_nof;

			//std::cout << info.hist_bef.sp << " " << info.hist_bef.step << " " << info.hist_bef.oc.size() << " " << info.hist_bef.oc[0] << " " << info.hist_bef.ocm << "\n";
			//std::cout << info.hist2_bef.sp << " " << info.hist2_bef.step << " " << info.hist2_bef.oc.size() << " " << info.hist2_bef.oc[0] << " " << info.hist2_bef.ocm << "\n";


			histogram hist_data_vol, hist_data_nof;
			hist_data_vol.read_hist_vol();
			hist_data_nof.read_hist_nof();

			std::cout << conp.nxyz << " \n";
			std::cout << "Removable: " << LAG_removable(conp, conpc, vec, hpar_estim, info, hist_data_nof, hist_data_vol, 0, 1, 0, 1, 0, 1) << " \n";

			std::cout << vec.size() << " " << vec.size()/3 << " \n";
			for (int ii = 0; ii < 10; ii++) { std::cout << vec[3*ii] << " " << vec[3 * ii + 1] << " " << vec[3 * ii + 2] << " \n"; } 	

			theta_e = 0, zet_e = 0;
			int N = 1000;
			LAG_estim(theta_e, zet_e, conp, conpc, N, hpar_estim, info, hist_data_nof, hist_data_vol, 0, 1, 0, 1, 0, 1);

			std::cin >> i;
			return 0;

			//std::cout << hist_data_nof.sp << " " << hist_data_nof.step << " " << hist_data_nof.oc.size() << " " << hist_data_nof.oc[0] << " " << hist_data_nof.ocm << "\n";
			//std::cout << hist_data_vol.sp << " " << hist_data_vol.step << " " << hist_data_vol.oc.size() << " " << hist_data_vol.oc[0] << " " << hist_data_vol.ocm << "\n";

			//double test1 = hist_disp(info.hist_bef, hist_data_nof);
			//double test2 = hist_disp(info.hist2_bef, hist_data_vol);
			//if (test1 == 0) { std::cout << "hist_disp = 0 \n"; }
			//else { std::cout << test1 << " \n"; }
			//if (test2 == 0) { std::cout << "hist_disp = 0 \n"; }
			//else { std::cout << test2 << " \n"; }

			std::cout << "proporcni rozdil: " << hist_disp(info.hist_bef, hist_data_nof) << " absolutni rozdil: " << hist_dis(info.hist_bef, hist_data_nof) << "\n";
			std::cout << "proporcni rozdil: " << hist_disp(info.hist2_bef, hist_data_vol) << " absolutni rozdil: " << hist_dis(info.hist2_bef, hist_data_vol) << "\n";
			


			std::cout << " Histogram volume: \n";
			std::cout << hist_vol.oc.size() << " " << hist_data_vol.oc.size() << " \n";
			std::cout << hist_vol.so << " " << hist_data_vol.so << " \n";
			std::cout << hist_vol.noc << " \n";
			for (int ii = 0; ii < hist_vol.oc.size(); ii++) { std::cout << hist_vol.oc[ii] << " "; } 	std::cout << "\n";
			for (int ii = 0; ii < hist_data_vol.oc.size(); ii++) { std::cout << hist_data_vol.oc[ii] << " "; } 	std::cout << "\n";

			std::cout << " Histogram nof: \n";
			std::cout << hist_nof.oc.size() << " " << hist_data_nof.oc.size() << " \n";
			std::cout << hist_nof.so << " " << hist_data_nof.so << " \n";
			std::cout << hist_nof.noc << " \n";
			for (int ii = 0; ii < hist_nof.oc.size(); ii++) { std::cout << hist_nof.oc[ii] << " "; } 	std::cout << "\n";
			for (int ii = 0; ii < hist_data_nof.oc.size(); ii++) { std::cout << hist_data_nof.oc[ii] << " "; } 	std::cout << "\n";

			bool cell;
			voronoicell c;
			double xn0 = 0, xn1 = 0;
			for (int j = 0; j < conp.nxyz; j++) { // loop over boxes
				for (int i = 0; i < conp.co[j]; i++) { // loop over particles in considered box

					cell = conp.compute_cell(c, j, i);
					if (cell == true) {
						xn0 = xn0 + c.number_of_faces();		// nof
						xn1 = xn1 + c.volume();					// volume
					} // the cell was computed 
				}
			}
			std::cout << " Prumer: " << (xn0 / info.tp_bef) << " ; " << info.mean_bef[0] / info.tp_bef << "; " << info.mean_bef[0] << " " << info.mean_aft[0] << "\n";
			std::cout << " Prumer: " << (xn1 / info.tp_bef) << " ; " << info.mean_bef[1] / info.tp_bef << "; " << info.mean_bef[1] << " " << info.mean_aft[1] << "\n";

			//FILE *f;
			//f = fopen("run.txt", "a"); // run.txt se otevre k zapisu na konec souboru

			std::cout << " BDMA \n";
            //h = 3000;
			h = 1000000;
            std::cout << " Number of steps: " << h << " \n";
			i = 0;

			double ar = 0;

			for (g = 0; g < h; g++) {
                                //std::cout << " STEP " << i++ << '\n';
				if(conp.total_particles() < 100){
					std::cout << "Vycerpani castic: " << n_add << " " << n_del << " " << n_mov << " " << n_part << "\n";
					std::cout << "step " << g-1 << "\n";
					std::cout << "conp.total_particles() = " << conp.total_particles() << " ; n = " << n << "\n";

					return 0;
				}
				// fci LAG_bdma se musi predat dve totozne kopie containeru (conp, conpc); fce provede zmeny v conpc za ucelem prepoctu energie, tyto zmeny
				//		je pak nutne vratit zpet a navic ztotoznit strukturu containeru, byla-li zmenena
				LAG_bdma(n, fid, conp, conpc, sigma, theta, zet, hpar, n_add, n_del, n_mov, info, hist_data_nof, hist_data_vol);
				//bdma_step(n, fid, conp, sigma, theta[0], zet, hpar[0], hpar[1], hpar[2], hpar[3], n_add, n_del, n_mov, ar);

				if (((g + 1) % 1000) == 0) {
					n_part = n_part + n_add - n_del; // aktualizace celkoveho poctu castic
					std::cout << "STEP " << g+1 << " : " << n_add << " " << n_del << " " << n_mov << " " << n_part << " (" << empty_cells(conp) << ") " << "\n";
					//std::cout << "     " << info.tp_bef << " ; " << info.mean_bef[0] / info.tp_bef << " ; " << info.mean_bef[1] / info.tp_bef << " \n";
					//std::cout << "     " << info.mean_bef[0] / info.tp_bef << " ; " << hist_disp(info.hist_bef, hist_nof) << "  " << hist_dis(info.hist_bef, hist_nof) << "\n";
					//std::cout << "     " << info.mean_bef[0] / info.tp_bef << " " << V_function(info.hist_bef, hist_data_nof) << " ; " << hist_disp(info.hist_bef, hist_data_nof) << "  " << hist_dis(info.hist_bef, hist_data_nof) << "\n";
					//std::cout << "     " << info.mean_bef[1] / info.tp_bef << " " << V_function_2(info.hist2_bef, hist_data_vol) << " ; " << hist_disp(info.hist2_bef, hist_data_vol) << "  " << hist_dis(info.hist2_bef, hist_data_vol) << "\n";

					//std::cout << "var " << info.var_bef[0] / (info.tp_bef-1) << " \n";
					//for (int ii = 0; ii < info.hist_bef.oc.size(); ii++) { std::cout << info.hist_bef.oc[ii] << " "; } 	std::cout << "\n";
					//std::cout << info.hist_bef.noc << " \n";
					n_add = 0;
					n_del = 0;
					n_mov = 0;

				}

			//	if (((g + 1) % 2000) == 0) {
			//	-	std::cout << hist_disp(info.hist_bef, hist_nof) << " \n";
			//		for (int ii = 0; ii < info.hist_bef.oc.size(); ii++) { std::cout << info.hist_bef.oc[ii] << " "; } 	std::cout << "\n";
			//		std::cout << info.hist_bef.noc << " \n";
			//	-	for (int ii = 0; ii < hist_nof.oc.size(); ii++) { std::cout << hist_nof.oc[ii] << " "; } 	std::cout << "\n";
			//	}

				//std::cout << n_add << " " << n_del << " " << n_mov << " " << n_chr << " " << n_part << "\n";

				//n_add = 0;
				//n_del = 0;
				//n_mov = 0;
				//n_chr = 0;
				  
			}
			std::cout << "END: " << n_add << " " << n_del << " " << n_mov << " " << conp.total_particles() << "\n";

			n_add = 0;
			n_del = 0;
			n_mov = 0;

		
			//if (((g+1) % 1000000) == 0) { // kazdy milionty krok se vypisi zakladni statistiky
			//	i = (g + 1) / 1000000;
			//	cell_stats(con,i);
			//	face_stats(con,i);
			//	std::cout << i << " \n";

			std::cout << "Empty cells: " << empty_cells(conp) << "\n";


			xn0 = 0; xn1 = 0;
			for (int j = 0; j < conp.nxyz; j++) { // loop over boxes
				for (int i = 0; i < conp.co[j]; i++) { // loop over particles in considered box

					cell = conp.compute_cell(c, j, i);
					if (cell == true) {
						xn0 = xn0 + c.number_of_faces();		// nof
						xn1 = xn1 + c.volume();					// volume
					} // the cell was computed 
				}
			}
			std::cout << " Prumer: " << (xn0 / info.tp_bef) << " ; " << info.mean_bef[0] / info.tp_bef << "; " << info.mean_bef[0] << " " << info.mean_aft[0] << "\n";
			std::cout << " Prumer: " << (xn1 / info.tp_bef) << " ; " << info.mean_bef[1] / info.tp_bef << "; " << info.mean_bef[1] << " " << info.mean_aft[1] << "\n";
			
			cell_stats(conp);
			face_stats(conp);

			write_container(conp);

			std::cout << info.tp_bef << "\n";

			std::cout << "proporcni rozdil: " << hist_disp(info.hist_bef, hist_data_nof) << " absolutni rozdil: " << hist_dis(info.hist_bef, hist_data_nof) << "\n";
			std::cout << "proporcni rozdil: " << hist_disp(info.hist2_bef, hist_data_vol) << " absolutni rozdil: " << hist_dis(info.hist2_bef, hist_data_vol) << "\n";


			if (feasibility(conp, hpar[0], hpar[1], hpar[2])) { std::cout << " Feasible \n"; }
			else { std::cout << " Unfeasible \n"; }

			hardcore_estim(conp, alfa_e, beta_e, B_e);
			std::cout << alfa_e << " " << beta_e << " " << B_e << " \n";

			std::cin >> i;
			
 
	return 0;
}
