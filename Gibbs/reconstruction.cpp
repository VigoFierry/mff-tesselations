#include "Header.h"

using namespace voro;

// reconstruction of the tessellation



// the intensity is given at the beginning, we need to create starting configuration


void reconstruct_naive() {

	int N = 1057; //intensity

	//create containers:
	double ax = 0, bx = 1, ay = 0, by = 0, az = 0, bz = 1; //rozmery kontejneru
	int nx, ny, nz; //pocty boxu

	// idealne 5 castic na box --> pocet boxu = N/5 
	// viz void pre_container_base::guess_optimal(int &nx,int &ny,int &nz) for inspiration
	double dx = bx - ax, dy = by - ay, dz = bz - az;
	double ilscale = pow(N / (optimal_particles*dx*dy*dz), 1 / 3.0);
	nx = int(dx*ilscale + 1);
	ny = int(dy*ilscale + 1);
	nz = int(dz*ilscale + 1);

	// Create a container with the geometry given above, and make it periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container_poly con(ax, bx, ay, by, az, bz, nx, ny, nz, true, true, true, 8);
	container_poly con_copy(ax, bx, ay, by, az, bz, nx, ny, nz, true, true, true, 8);

	double x, y, z, r;
	double MR = 0.2; //max radius
	for (int i = 0; i<N; i++) { // Randomly add particles into the container
		x = ax + uniform(0, 1)*(bx - ax);
		y = ay + uniform(0, 1)*(by- ay);
		z = az + uniform(0, 1)*(bz - az);
		r = uniform(0, 1)*MR;
		con.put(i, x, y, z, r);
		con_copy.put(i, x, y, z, r);
	}

	histogram hist_data_vol, hist_data_nof;
	hist_data_vol.read_hist_vol();	//nacteni histogramu
	hist_data_nof.read_hist_nof();

	con_info info;
	info.get_tp(con);
	info.get_meansum(con);

	histogram hist_vol, hist_nof;
	hist_nof.create_hist_int(con); // nof
	hist_vol.create_hist_double(con); // volume
									   //std::cout << hist.sp << " " << hist.step << " " << hist.oc.size() << " " << hist.oc[0] << " " << hist.ocm << "\n";
	info.hist2_bef = hist_vol;
	info.hist2_aft = hist_vol;
	info.hist_bef = hist_nof;
	info.hist_aft = hist_nof;


	//navrh zmeny:
	double xx, yy, zz, rr;	//nove hodnoty
	int ijk, q;
	int id = uniform_int(1,N);
	find_pos(ijk, q, id, &con);					
	x = con.p[ijk][4 * q]; y = con.p[ijk][4 * q + 1]; z = con.p[ijk][4 * q + 2]; r = con.p[ijk][4 * q + 3];
	xx = ax + uniform(0, 1)*(bx - ax);
	yy = ay + uniform(0, 1)*(by - ay);
	zz = az + uniform(0, 1)*(bz - az);
	rr = uniform(0, 1)*MR;

	erase(ijk, q, &con_copy);				// smazani castice
	con_copy.put(id, xx, yy, zz, rr);		// pridani nove castice (ID zustava zachovano)

	LAG_container(con, con_copy, 3, id);
	std::vector<int> cells; std::vector<int> cells_pos;
	LAG_cells(con, con_copy, 3, id, info, cells, cells_pos);
	std::vector<int> sec; std::vector<int> sec_pos;
	LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);
	bool pripustnost = true;
	//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos);

	long double val_bef, val_aft;
	int i, ijk_mb, q_mb;
	bool cell, cell2;
	voronoicell_neighbor c;
	long double val1_a = 0, val1_b = 0, val2_a = 0, val2_b = 0;
	int tp = 0;
	find_pos(ijk_mb, q_mb, id, &con);

	bool change = 1;
	if (pripustnost == false) { change = 0; }
	else {

		// stejny cyklus jako pro zjisteni secondary neighbours
		//# pragma omp parallel for shared(con, newcon)
		for (i = 0; i < cells.size(); i++) {
			val_bef = 0; val_aft = 0;
			if (cells[i] == id) {
				
//				if (type == 3) {		// move
					cell = con.compute_cell(c, ijk_mb, q_mb);
					if (cell == true) {
						tp--;
						// 1. volume
						val_bef = c.volume();
						//val2_b = val2_b + val_bef;
						info.hist2_aft.hist_act(val_bef, 0);
						//					part2 = part2 + V_function_2(val_bef, hist2);
						// 2. number of faces
						val_bef = c.number_of_faces();
						val1_b = val1_b + val_bef;
						info.hist_aft.hist_act(val_bef, 0);
						//					part1 = part1 + V_function(val_bef, hist);

					}
					cell = con_copy.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
					if (cell == true) {
						tp++;
						// 1. volume
						val_aft = c.volume();
						//val2_a = val2_a + val_aft;
						info.hist2_aft.hist_act(val_aft, 1);
						//					part2 = part2 - V_function_2(val_aft, hist2);
						// 2. number of faces
						val_aft = c.number_of_faces();
						val1_a = val1_a + val_aft;
						info.hist_aft.hist_act(val_aft, 1);
						//					part1 = part1 - V_function(val_aft, hist);
					}

//				}
			}	// end if (changed particle)					// zde nutne rozlisit pripady add/delete/move
			else {
				cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					tp--;
					// 1. volume	
					val_bef = c.volume();
					//val2_b = val2_b + val_bef;
					info.hist2_aft.hist_act(val_bef, 0);
					//				part2 = part2 + V_function_2(val_bef, hist2);
					// 2. number of faces
					val_bef = c.number_of_faces();
					val1_b = val1_b + val_bef;
					info.hist_aft.hist_act(val_bef, 0);
					//				part1 = part1 + V_function(val_bef, hist);

					// position of centroid
					// number of edges
					// surface area
					// total edge distance
					// max radius squared
					// face areas
					// face perimeters
					// face vertices
					// ...
					// 6. radius
					//val_bef = con.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];

				}
				cell2 = con_copy.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell2 == true) {
					tp++;
					// 1. volume	
					val_aft = c.volume();
					//val2_a = val2_a + val_aft;
					info.hist2_aft.hist_act(val_aft, 1);
					//				part2 = part2 - V_function_2(val_aft, hist2);
					// 2. number of faces
					val_aft = c.number_of_faces();
					val1_a = val1_a + val_aft;
					info.hist_aft.hist_act(val_aft, 1);
					//				part1 = part1 - V_function(val_aft, hist);

					// ...
					// 6. radius
					//val_aft = newcon.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];
				}
				//			std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm)  << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " " << val_aft << " \n"; ////////////
				//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
				//part1 = part1 + V_function(val_bef, hist) - V_function(val_aft, hist);
				// energy = energy + val_bef - val_aft;
				// energy = energy + min(val_bef, K) - min(val_aft, K);
				// energy = energy + (abs(val_bef - n0))/n0 - (abs(val_aft - n0))/n0;

				// energy = energy + V_function(val_bef) - V_function(val_aft);
			} // end else (changed particle)

			  //std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////
			  //energy = energy + (1 - hist.hist_value(val_bef) / hist.ocm) - (1 - hist.hist_value(val_aft) / hist.ocm);	// ____________________ znamenka
		}
	}






	if (change) {																// ZMENA PROVEDENA
		erase(ijk, q, &con);														// smazani castice
		con.put(id, xx, yy, zz, rr);												// pridani nove castice (ID zustava zachovano)
		info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
																//info.mean_bef[1] = info.mean_aft[1]; 
																//info.var_bef[0] = info.var_aft[0];
		info.tp_bef = info.tp_aft;
		info.hist_bef = info.hist_aft;
		info.hist2_bef = info.hist2_aft;
		//			std::cout << "YES \n";
//		rn_mov++;																	// counter of moved particles
	}
	else {																			// ZMENA ZAMITNUTA - vrat conp_copy do puvodniho stavu
		find_pos(ijk, q, id, &con_copy);
		erase(ijk, q, &con_copy);													// smazani castice
		con_copy.put(id, x, y, z, r);												// pridani nove castice (ID zustava zachovano)

		info.hist_aft = info.hist_bef;
		info.hist2_aft = info.hist2_bef;
		info.tp_aft = info.tp_bef;
		info.mean_aft[0] = info.mean_bef[0];
		//info.mean_aft[1] = info.mean_bef[1];
		//			std::cout << "NO \n";
	}
}