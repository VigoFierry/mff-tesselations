#include "Header.h"


using namespace voro;


void LAG_bdma_cout(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, std::vector<double> theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info, histogram &hist, histogram &hist2)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in]			con							the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in]			conp_copy					copy of the original container (as reference - will be manipulated).
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						vector of smoothing parameters
	// [in]			zet							intensity constant of reference Poisson process
	// [in]			hard_par					(= alpha, beta, B, iota) vector ofhardcore parameters
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles

	// schema:
	// 1) navrhni zmenu - add/delete/move
	// 2) rozdvoj container a v kopii proved zmenu
	// 3) predej informace (dva containery, typ operace, id, parametry) fci LAG_recompute
	// 4) na zaklade spoctene energie rozhodni o provedeni operace
	// 5) aktualizuj container 

	double rn1 = uniform(0, 1);   // random number between 0 and 1  - navrh zmeny
	double rn2 = uniform(0, 1);   // random number between 0 and 1 

								  //rn2 = 0.0000000000000001;
	int typ; // 1 - add, 2 - delete, 3 - move

			 //std::cout << "jsem tu \n";

			 // kopie containeru:
			 // container_poly conp_copy; ................................. no default constructor exists !!!
			 //container_poly conp_copy(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8); 
			 //conp_copy = conp;
			 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			 // - vyreseno predanim dvou identickych containeru teto funkci, predan je original referenci - nebude s nim manipulovano, a kopie opet referenci - poslouzi 
			 //		k vypoctum (predavani jako lokalni promenna ma za nasledek spadnuti programu)


	double nx, ny, nz, x, y, z, r, nr;
	double pst;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	//rn1 = 0.2;
	if (rn1 <= (const1)) {
		typ = 1;
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// the radius of the new particle can be generated from the prespecified distribution:
		r = 0.05*uniform(0, 1);
		//r = gamma(3, 0.5);
		//r = triangle(0.005, 0.03, 0.02);
		// or using the average radius in container:
		//r = ave_rad(conp);

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
																	// no = empty_cells(conp);

		conp_copy.put(id, nx, ny, nz, r);							// pridani castice do kopie

		//		pst = LAG_recompute(conp, conp_copy, typ, id, hard_par, theta, info);	// prepocet energie
		//		pst = try_add(id, nx, ny, nz, r, conp, theta, alfa, beta, B, iota);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, id);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos);


		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos); 
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			LAG_V1(conp, conp_copy, typ, id, info, hist, hist2, parts, cells, cells_pos); // hist se p
			pst = exp(theta[0] * parts[0]); //+ theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, id, info, cells, cells_pos, sec, sec_pos);
			pst = pst*exp(theta[2] * parts[2]);
			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
		}

		// pst = 0.5;
		std::cout << pst << " ";
		pst = pst*(zet / (no + 1));
// -		pst = pst / (no + 1); // without z
		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {						// ZMENA PROVEDENA - proved zmenu i v con
			conp.put(id, nx, ny, nz, r);		// pridani castice
												// can NOT be simply overwritted .......... conp = &conp_copy;

			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
																	//info.mean_bef[1] = info.mean_aft[1];
																	//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
			std::cout << "YES \n";
			if (id > npart) {					// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;				// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {								// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;							// counter of added particles

												//maintaince of empty cells
												//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
												//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
												//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
												//			}
		}
		else {									// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			find_pos(ijk, q, id, &conp_copy);
			erase(ijk, q, &conp_copy);

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
			std::cout << "NO \n";
		}

//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		typ = 2;
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; r = conp.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

		no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
																	// no = empty_cells(conp);

		erase(ijk, q, &conp_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);	// prepocet energie
//		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);	// spocte pravdepodobnost se kterou dojde k operaci DELETE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos);
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			LAG_V1(conp, conp_copy, typ, del, info, hist, hist2, parts, cells, cells_pos); // hist se p
			pst = exp(theta[0] * parts[0]); // +theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, del, info, cells, cells_pos, sec, sec_pos);
			pst = pst*exp(theta[2] * parts[2]);
			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
		}


		// pst = 0.5;
		std::cout << pst << " ";
		pst = pst*(no / zet);
		// -		pst = pst * no; // without z
		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {											// ZMENA PROVEDENA - proved zmenu i v con
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
																	//info.mean_bef[1] = info.mean_aft[1];
																	//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
			std::cout << "YES \n";
			rn_del++;												// counter of deleted particles

																	//maintaince of empty cells
																	//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
																	//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
																	//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
																	//			}
		}
		else {														// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			conp_copy.put(del, x, y, z, r);

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
			std::cout << "NO \n";
		}
	}     // DEATH

	if (((const2) < rn1)) {
		typ = 3;
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
									// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
		// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = 0.05*uniform(0, 1);
		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		erase(ijk, q, &conp_copy);														// smazani castice
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);						// prepocet energie
//		pst = try_MOVE(ijk, q, nx, ny, nz, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos);
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			LAG_V1(conp, conp_copy, typ, del, info, hist, hist2, parts, cells, cells_pos); // hist se p
			pst = exp(theta[0] * parts[0]); // +theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, del, info, cells, cells_pos, sec, sec_pos);
			pst = pst*exp(theta[2] * parts[2]);
			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
		}


		// pst = 0.5;
		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {																// ZMENA PROVEDENA
			erase(ijk, q, &conp);														// smazani castice
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
																	//info.mean_bef[1] = info.mean_aft[1]; 
																	//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
			std::cout << "YES \n";
			rn_mov++;																	// counter of moved particles

																						//maintaince of empty cells
																						//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
																						//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
																						//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
																						//			}
		}
		else {																			// ZMENA ZAMITNUTA - vrat conp_copy do puvodniho stavu
			find_pos(ijk, q, del, &conp_copy);
			erase(ijk, q, &conp_copy);													// smazani castice
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
			std::cout << "NO \n";
		}

	}     // MOVE 

}

void LAG_bdma(long &npart, std::vector<int> &fid, voro::container_poly &conp, voro::container_poly &conp_copy, double sigma, std::vector<double> theta, double zet, std::vector<double> hard_par, long &rn_add, long &rn_del, long &rn_mov, con_info &info, histogram &hist, histogram &hist2)
{
	// [in,out]		npart						number of particles in the container.
	// [in,out]		fid							vector of available id´s.
	// [in]			con							the original container with stored particles (as reference - will not be changed, but the structure can be rearranged).
	// [in]			conp_copy					copy of the original container (as reference - will be manipulated).
	// [in]			sigma						parameter of normal distribution.
	// [in]			theta						vector of smoothing parameters
	// [in]			zet							intensity constant of reference Poisson process
	// [in]			hard_par					(= alpha, beta, B, iota) vector ofhardcore parameters
	// [out]		n_add, n_del, n_mov, n_chr	counters of added/deleted/moved/radius changed particles

	// schema:
	// 1) navrhni zmenu - add/delete/move
	// 2) rozdvoj container a v kopii proved zmenu
	// 3) predej informace (dva containery, typ operace, id, parametry) fci LAG_recompute
	// 4) na zaklade spoctene energie rozhodni o provedeni operace
	// 5) aktualizuj container 

	double rn1 = uniform(0, 1);   // random number between 0 and 1  - navrh zmeny
	double rn2 = uniform(0, 1);   // random number between 0 and 1 

								  //rn2 = 0.0000000000000001;
	int typ; // 1 - add, 2 - delete, 3 - move

			 //std::cout << "jsem tu \n";

			 // kopie containeru:
			 // container_poly conp_copy; ................................. no default constructor exists !!!
			 //container_poly conp_copy(0, 1, 0, 1, 0, 1, nx, ny, nz, true, true, true, 8); 
			 //conp_copy = conp;
			 // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			 // - vyreseno predanim dvou identickych containeru teto funkci, predan je original referenci - nebude s nim manipulovano, a kopie opet referenci - poslouzi 
			 //		k vypoctum (predavani jako lokalni promenna ma za nasledek spadnuti programu)


	double nx, ny, nz, x, y, z, r, nr;
	double pst = 1;
	double const1 = 1;
	double const2 = 2;
	double const3 = 3;

	int del, id, no;
	int ijk, q; // ijk_add, q_add; unsigned int i,j,k; 

				//std::cout << "pst: " << rn1 << " " << rn2 << "\n";
	const1 = const1 / const3;
	const2 = const2 / const3;

	// std::cout << "fid: ";	///////////////////////////////////////////////////////////////////////////////////////////
	// for (i = 0; i < fid.size(); i++) {  ///////////////////////////////////////////////////////////////////////////////
	// std::cout << fid[i] << " ";
	// } std::cout << " \n";

	//rn1 = 0.8;
	if (rn1 <= (const1)) {
		typ = 1;
		nx = uniform(0, 1); ny = uniform(0, 1); nz = uniform(0, 1); // coordinates of new particle
																	// the radius of the new particle can be generated from the prespecified distribution:
		r = 0.05*uniform(0, 1);
		//r = gamma(3, 0.5);
		//r = triangle(0.005, 0.03, 0.02);
		// or using the average radius in container:
		//r = ave_rad(conp);

		if (fid.size() > 0) {										// najiti vhodneho ID
			id = fid[0];
		}   // zamerem bylo zjistit zdali je ve fid ulozene nejake volne id ci ne???
		else { id = npart + 1; }  // pokud nebylo dostupne ID ve fid, pouzije se maximalni hodnota ID + 1
//		std::cout << "ADD " << id << " " << nx << " " << ny << " " << nz << " " << r << " \n"; // popis cinosti

		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky
																	// no = empty_cells(conp);
		no = info.tp_bef;

		conp_copy.put(id, nx, ny, nz, r);							// pridani castice do kopie

//		pst = LAG_recompute(conp, conp_copy, typ, id, hard_par, theta, info);	// prepocet energie
//		pst = try_add(id, nx, ny, nz, r, conp, theta, alfa, beta, B, iota);		  	// spocte pravdepodobnost se kterou dojde k operaci ADD
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, id);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, id, info, cells, cells_pos);

		 
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, id, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos); 
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			//LAG_V1(conp, conp_copy, typ, id, info, hist, hist2, parts, cells, cells_pos); // hist se p
			//pst = exp(theta[0] * parts[0] +theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, id, info, cells, cells_pos, sec, sec_pos);
			pst = exp(theta[2] * parts[2]);
			//pst = pst*exp(theta[2] * parts[2]);
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
			//std::cout << "ADD " << parts[2] << " " << pst << "\n";
		} 

	 // pst = 0.5;
//		std::cout << pst << " ";
	//	std::cout << V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist) << "\n";
		pst = pst*(zet / (no + 1)); 
// -		pst = pst / (no + 1); // without z
//		std::cout << pst << " ";

		// rn = rnd();							// generate random number between 0 and 1 
		if (rn2 < pst) {						// ZMENA PROVEDENA - proved zmenu i v con
			conp.put(id, nx, ny, nz, r);		// pridani castice
												// can NOT be simply overwritted .......... conp = &conp_copy;
			
			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
			//info.mean_bef[1] = info.mean_aft[1];
			//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
//			std::cout << "YES \n";
			if (id > npart) {					// pouzilo-li se ID o jedna vyssi nez maximalni ID v kontejneru a pridani
				npart = npart + 1;				// castice se uskutecni, pak je nutne navysit toto maximalni ID o 1
			}
			else {								// pouzilo-li se ID z fid, musi ve fid smazat
				fid.erase(fid.begin());
			}
			rn_add++;							// counter of added particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//			}
		}
		else {									// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			find_pos(ijk, q, id, &conp_copy);
			erase(ijk, q, &conp_copy);

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
//			std::cout << "NO \n";
		}

		//		std::cout << npart << "\n";  /////////////////////////////////////////////////////////////////////////////////
	}     // BIRTH

	if (((const1) < rn1) && (rn1 <= (const2))) {
		typ = 2;
		del = uniform_int(1, conp.total_particles());	// vybere castici (resp. jeji poradi v containeru) ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; r = conp.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

		//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
																	// no = empty_cells(conp);
		no = info.tp_bef;

		erase(ijk, q, &conp_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);	// prepocet energie
		//		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);	// spocte pravdepodobnost se kterou dojde k operaci DELETE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos);
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			//LAG_V1(conp, conp_copy, typ, del, info, hist, hist2, parts, cells, cells_pos); // hist se p
			//pst = exp(theta[0] * parts[0] +theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, del, info, cells, cells_pos, sec, sec_pos);
			pst = exp(theta[2] * parts[2]);
			//pst = pst*exp(theta[2] * parts[2]);
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
			//std::cout << "DELETE " << parts[2] << " " << pst << "\n";
		}

		
		// pst = 0.5;
//		std::cout << pst << " ";
		pst = pst*(no / zet); 
// -		pst = pst * no; // without z
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {											// ZMENA PROVEDENA - proved zmenu i v con
			erase(ijk, q, &fid, &conp);								// smazani castice - uvolni ID do fid k opetovnemu pouziti
			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
			//info.mean_bef[1] = info.mean_aft[1];
			//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
//			std::cout << "YES \n";
			rn_del++;												// counter of deleted particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//			}
		}
		else {														// ZMENA ZAMITNUTA - navrat con_copy do puvodniho stavu
			conp_copy.put(del, x, y, z, r);

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
//			std::cout << "NO \n";
		}
	}     // DEATH

	if (((const2) < rn1)) {
		typ = 3;
		del = uniform_int(1, conp.total_particles());	// vybere castici ke smazani (uniforme z celkoveho poctu castic)
														// del neni ID!!!
		find_part(ijk, q, del, &conp);					// najde tuto castici, tj jeji polohu ve strukture containeru 
		del = conp.id[ijk][q];							// do del se nyni ulozi ID teto castice misto jejiho poradi
		x = conp.p[ijk][4 * q]; y = conp.p[ijk][4 * q + 1]; z = conp.p[ijk][4 * q + 2]; // urci jeji souradnice
		r = conp.p[ijk][4 * q + 3];	// a jeji radius
									// vygeneruje novou castici z normalniho rozdeleni:
		nx = normal(x, sigma); ny = normal(y, sigma); nz = normal(z, sigma);			// coordinates of new particle
//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti
																															// these coordinates could be outside the window [0,1]^3 and in this case it is necessary to remap it into window:
		nx = nx - step_int(nx); ny = ny - step_int(ny); nz = nz - step_int(nz);
		// novy polomer nezavisi na to puvodnim !!! tj z hlediska polomeru nejde o "move"
		nr = 0.05*uniform(0, 1);
		//		std::cout << "MOVE " << del << " " << x << " " << y << " " << z << " --> " << nx << " " << ny << " " << nz << "\n"; // popis cinosti

		erase(ijk, q, &conp_copy);														// smazani castice
		conp_copy.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)

//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);						// prepocet energie
		//		pst = try_MOVE(ijk, q, nx, ny, nz, nr, conp, theta, alfa, beta, B, iota);		// urci pravdepodobnost se kterou dojde k operaci MOVE
		// fci LAG_recompute si lze "poskladat":
		LAG_container(conp, conp_copy, typ, del);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(conp, conp_copy, typ, del, info, cells, cells_pos);
		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(conp, conp_copy, del, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos);
		if (pripustnost == false) { pst = 0; }
		else {
			std::vector<double> parts; parts.resize(3);
			//LAG_V1(conp, conp_copy, typ, del, info, hist, hist2, parts, cells, cells_pos); // hist se p
			//pst = exp(theta[0] * parts[0] +theta[1] * parts[1]);
			parts[2] = LAG_V2(conp, conp_copy, typ, del, info, cells, cells_pos, sec, sec_pos);
			pst = exp(theta[2] * parts[2]);
			//pst = pst*exp(theta[2] * parts[2]);
//			std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";
			//std::cout << "MOVE " << parts[2] << " " << pst << "\n";
		}

		
		// pst = 0.5;
//		std::cout << pst << " ";

		// rn = rnd();  // generate random number between 0 and 1 
		if (rn2 <= pst) {																// ZMENA PROVEDENA
			erase(ijk, q, &conp);														// smazani castice
			conp.put(del, nx, ny, nz, nr);												// pridani nove castice (ID zustava zachovano)
			info.mean_bef[0] = info.mean_aft[0];					// aktualizace infa
			//info.mean_bef[1] = info.mean_aft[1]; 
			//info.var_bef[0] = info.var_aft[0];
			info.tp_bef = info.tp_aft;
			info.hist_bef = info.hist_aft;
			info.hist2_bef = info.hist2_aft;
//			std::cout << "YES \n";
			rn_mov++;																	// counter of moved particles

			//maintaince of empty cells
//			for (int iu = 0; iu < info.empty.size() / 2; iu++) {
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &fid, &conp);
//				erase(info.empty[2 * iu], info.empty[2 * iu + 1], &conp_copy);
//			}
		}
		else {																			// ZMENA ZAMITNUTA - vrat conp_copy do puvodniho stavu
			find_pos(ijk, q, del, &conp_copy);
			erase(ijk, q, &conp_copy);													// smazani castice
			conp_copy.put(del, x, y, z, r);												// pridani nove castice (ID zustava zachovano)

			info.hist_aft = info.hist_bef;
			info.hist2_aft = info.hist2_bef;
			info.tp_aft = info.tp_bef;
			info.mean_aft[0] = info.mean_bef[0];
			//info.mean_aft[1] = info.mean_bef[1];
//			std::cout << "NO \n";
		}

	}     // MOVE 

}

// specification of the potential (outsourced energy)
double V_function(double val)
{
	//	[in]	val		double / int

	// volume					double
	// number of faces			int
	// position of centroid		double, double, double
	// number of edges			int
	// surface area				double
	// total edge distance		double
	// max radius squared		double
	// face areas				vector<double>
	// face perimeters			vector<double>
	// face vertices			vector<int>
	// ...
	// radius					double

	double K;
	K = 1;
	int n0;
	n0 = 14;

	// energy = energy + val_bef - val_aft;
	// energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	
	// energy = energy + min(val_bef, K) - min(val_aft, K);
	// energy = energy + (abs(val_bef - n0))/n0 - (abs(val_aft - n0))/n0;

	return val;
	min_max(val, K); return val;
	//return (1 - hist.hist_value(val) / hist.ocm);
	//return sqrt(1 - hist.hist_value(val) / hist.ocm);
	//return sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm));
	//return sqrt(sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm)));
	return ((abs(val - n0)) / n0);
	return sqrt((abs(val - n0)) / n0);

}

double V_function(double val, histogram &hist)
{
	//	[in]	val		double / int

	
	//double K = 1;
	double n0 = 12;
	//double n1 = 2;
	double ns2 = 100000;
	//double c0 = 1057;
	//double c = 2000;
	//double r0 = 0.005;

	// energy = energy + val_bef - val_aft; 
	//return val; 
	//min_max(val, K); return val; 
	//return ((pow(2 * ((1 - hist.hist_value(val) / hist.ocm) - 0.5), 5) / 2) + 0.5)*c;
	//return pow((1 - hist.hist_value(val) / hist.ocm),4);
	//return (1 - hist.hist_value(val) / hist.ocm); 
	//return sqrt((1 - hist.hist_value(val) / hist.ocm )); 
	//return sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm)); 
//	std::cout << hist.hist_value(val) << " " << hist.ocm << " " << hist.hist_value(val) / hist.ocm << " \n"; 
	//return sqrt(sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm))); 
	//return ((abs(val - n0)) / n0); 
	//return sqrt((abs(val - n0)) / n1); 
	//return abs(val - c0) / c;	// ... experiment pro souhrnnou energiu (total particles)
	// val = MEAN (vyberovy prumer)
	return sqrt((abs(val - n0))) * ns2;   // ... experiment pro souhrnnou energii (mean of nof)
	//return sqrt((abs(val - r0)) * ns2);   // ... experiment pro souhrnnou energii (var of volume)
	// val = vyberovy rozptyl

}

double V_function_2(double val, histogram &hist)
{
	//	[in]	val		double / int

	//double K = 1;
	double n0 = 0.0015;
	//double n1 = 2; 
	//double ns2 = 89442.7191;
	double ns2 = 10000000;
	//double c0 = 1057;
	//double c = 2000; 
	//double r0 = 0.005;

	// energy = energy + val_bef - val_aft; 
	//return val; 
	//min_max(val, K); return val; 
	//return ((pow(2 * ((1 - hist.hist_value(val) / hist.ocm) - 0.5), 5) / 2) + 0.5)*c;
	//return pow((1 - hist.hist_value(val) / hist.ocm),4);
	//return (1 - hist.hist_value(val) / hist.ocm); 
	//return sqrt((1 - hist.hist_value(val) / hist.ocm));
	//return sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm)); 
	//	std::cout << hist.hist_value(val) << " " << hist.ocm << " " << hist.hist_value(val) / hist.ocm << " \n"; 
	//return sqrt(sqrt(sqrt(1 - hist.hist_value(val) / hist.ocm))); 
	//return ((abs(val - n0)) / n0); 
	//return sqrt((abs(val - n0)) / n1); 
	//return abs(val - c0) / c;	// ... experiment pro souhrnnou energiu (total particles)
	// val = MEAN (vyberovy prumer)
	return sqrt((abs(val - n0))) * ns2;   // ... experiment pro souhrnnou energii (mean of vol)
	//return sqrt((abs(val - r0)) * ns2);   // ... experiment pro souhrnnou energii (var of volume)
	// val = vyberovy rozptyl

}

double V_function(histogram &val, histogram &hist)
{
	double k = 1000000; 
	double l = 1000000;
	 
	double dsc = hist_disp(val, hist);
	//return hist_dis(val, hist) * l; 
	// hist_dis vraci absolutni rozdil mezi histogramy - tj. krome tvaru histogramu je regulovan i pocet bunek (cetnosti)
	//return sqrt(hist_disp(val, hist)) * k;
	return sqrt(dsc) * pow(dsc * 100, 5);
	// hist_disp vraci relativni rozdil mezi histogramy - tj. zohlednuje se jen tvar, cetnosti jsou normovany (pomoci so nebo ocm)
}


double V_function_2(histogram &val, histogram &hist)
{
	double k = 1000000;
	double l = 1000000000;

	double dsc = hist_disp(val, hist);
	//return hist_dis(val, hist) * l; 
	// hist_dis vraci absolutni rozdil mezi histogramy - tj. krome tvaru histogramu je regulovan i pocet bunek (cetnosti)
	return sqrt(dsc) * pow(dsc * 100,5);
	// hist_disp vraci relativni rozdil mezi histogramy - tj. zohlednuje se jen tvar, cetnosti jsou normovany (pomoci so nebo ocm)
}

// pair potential
double V_function(double val1, double val2)
{
	return 0;
}
double V_function(double val1, double val2, histogram &hist)
{
	return 0;
}


double V2_function(voronoicell_neighbor &rc1, voronoicell_neighbor &rc2)
{
	// [in]		c1,c2				the cells for which is function V2 computed.
	
	// assumption: c1 and c2 are neighbouring cells

	// assumption: potential is bounded
	double K = 10;

	double a = rc1.volume();  // vypocet objemu
	double b = rc2.volume();
	//std::cout << "Volumes: " << a << " " << b << "\n";

	min_max(a, b);		// a >= b

	a = sqrt(a / b - 1);
	//std::cout << a << " ";
	min_max(a, K);
	
//	std::cout << " " << K ;
	return K;

}


double LAG_recompute(container_poly &con, container_poly &newcon, int type, int id, std::vector<double> &h_par, std::vector<double> &theta, con_info &info)
{
	//	[in]	con		container before the change
	//	[in]	newcon	container after the change (local variable, 
	//	[in]	type	suggested change (add/delete/move)
	//  [in]	id		change proposal (id of changed particle)
	//	[in]	h_par	vector of hardcore parameters (alpha, beta, B, iota, ...), typical length is 3 or 4
	//	[in]	theta	vector of smooth parameters (theta1, theta2, ...)

	// con je puvodni container

	// newcon je jeho kopie, v te provedeme zmenu (add/delete/move) - rekneme, ze zmena je jiz provedena
	// type reprezentuje zmenu: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)
	double energy = 0;


	// A) STRUCTURE OF CONS, structure must be kept identical (for both con and con_copy)
		// ADD - pridana castice je na konci nejakeho boxu v newcon, pozice ostatnich se shoduji
		// DELETE - castice se odstrani, pozice dalsich castic v boxu v newcon se posunou -> reseni: trik s odebranim a opetovnym pridanim
	if (type == 2) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}
	// MOVE - zmenena castice je v con i v newcon opet na ruzne pozici (move=delete+add), 
	// id castice bylo zachovano, ale pozice se zmeni, tim dojde k posunu pozic i ostatnich castic -> reseni: trik odeber a pridej
	if (type == 3) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}


// 0) ENERGY INDEPENDENT OF GEOMETRY - slepa vetev, udelano lip nize
	// - pokud je tento vypocet aktivni, pak je potreba v LAG_bdma zakomentovat casti s vyskytem parametru z
/*	double part0 = 0;
	int n0 = 1057;
	double cst = 1;
	cst = cst / 200;
	int sgnop=0;
	int tp = con.total_particles();  // !!!! total particles neq non-empty cells !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//std::cout << "tp: " << tp << " ";
	if (type == 1) { sgnop = -1; }
	if (type == 2) { sgnop = 1; }
	//std::cout << "sgnop: " << sgnop << " ";
	int sgn = 0;
	if (tp > n0) { sgn = 1; }
	if (tp < n0) { sgn = -1; }
	//std::cout << "sgn: " << sgn << " ";
	part0 = part0 + cst * abs(tp - n0) * sgn * sgnop;
	//std::cout << "comp: " << cst << " ";
	//	std::cout << "energy (intensity): " << part0 << " ";
	
	
	energy = theta[0] * part0;
	//std::cout << "expenergy: " << exp(energy) << " \n";
	return exp(energy);
*/


	// B) LIST OF MODIFIED PARTICLES
		// list ovlivnenych bunek se muze lisit pro oba kontejnery (pripady add a delete) - poznacit si odlisnost? - odlisnost v jedne castici (te pridavane/odebirane)
		// zaklad listu lisicich se bunek tvori sousede zmenene castice
	std::vector<int> cells;
	std::vector<int> vert;
	int ijk, q;
	int fng;
	int ijk_ma, q_ma, ijk_mb, q_mb; // a - after (newcon), b - before (con)
	voronoicell_neighbor c;
	bool cell;
	double x, y, z, xx, yy, zz;
	double x_a, y_a, z_a, x_b, y_b, z_b;
	//	double wx, wy, wz;
	bool in;
	int count;

	//std::vector<int> empty; // empty cells
	info.empty.clear();


// 1a) INITIAL PARTICLES
	if (type == 1) {
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell == true) {						//  muze se stat, ze pridany generator vytvori prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);
			
		} // otherwise mozaika se nezmenila = energie se nezmenila = return exp(0)
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
			return 1;
		}
	}

	if (type == 2) {
		find_pos(ijk_mb, q_mb, id, &con);
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze odebirany generator tvoril prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);

			// assumption: no empty cells in the previous configuration
		} // otherwise mozaika se nezmenila = energie se nezmenila = return 0
	}

	std::vector<int> neigh;
	bool cell2;

	if (type == 3) {
		find_pos(ijk_mb, q_mb, id, &con);				// pozice generatoru NEMUSI byt stejna v con i v newcon !!! (zmena souradnic muze zpusobit zmenu boxu)
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze generator tvoril prazdnou bunku !!! (assumption: there were no empty cells before)
			c.neighbors(cells);

		}
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell2 = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell2 == true) {					//  muze se stat, ze zmeneny generator vytvori prazdnou bunku !!!
			c.neighbors(neigh);

		}
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
		}

		merge(cells, neigh);

	}


	//std::cout << "jsem tu \n";
//	std::cout << cell << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	if (type == 3) { std::cout << cell2; } ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	for (int ii = 0; ii < cells.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << cells[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	for (int ii = 0; ii < io.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << io[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	for (int ii = 0; ii < ap.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << ap[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 1b) ALGORITHM for SEARCHING POSSIBLY MODIFIED PARTICLES
	// implementace algoritmu vyhledani zmenenych bunek
	// cells obsahuje ids bunek, ktere se nachazi jak v con tak i v newcon
	// funkce ktera porovna dva listy integeru a vrati jejich odlisnosti
	int i;
	std::vector<int> neighi, verti;
	std::vector<int> cells_pos;
	std::vector<bool> newi;
	bool shoda = false;
	// kontrolni mechanismus:
	int a, b;

	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			if (type == 2) { cells_pos.push_back(-1); cells_pos.push_back(-1); } // delete (castice neni v newcon)
			else { cells_pos.push_back(ijk_ma); cells_pos.push_back(q_ma); } // uloz polohu z newcon
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede uz byli uvazovani
		else {
			find_pos(ijk, q, cells[i], &con);		// pozice generatoru je stejna v con i v newcon 
			//std::cout << "con: " << ijk << " " << q << " "; //////////////////////////////////////////////////////////////////////////////////////////////////
			cells_pos.push_back(ijk); cells_pos.push_back(q);
			cell = con.compute_cell(c, ijk, q);
			if (cell == true) {						// ______________________ sousede bunky id v con, pokud existuji !!! 
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, ijk, q);
			if (cell2 == true) {					// ______________________ sousede bunky id v newcon, pokud existuji !!! 
				c.neighbors(neighi);

			}
			else {
				info.empty.push_back(ijk);
				info.empty.push_back(q);
			}

			// kontrola:
			find_pos(a, b, cells[i], &newcon);		// pozice generatoru je stejna v con i v newcon ////////////////////////////////////////////////////////////
			//std::cout << "& newcon: " << a << " " << b << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////
			if (a == ijk && b == q) {}
			else { std::cout << "ERROR: (LAG_recompute) not corresponding placement! \n"; } //////////////////////////////////////////////////////////


		} // END if..else (id)
		// najdi v cem se lisi neigh a neighi --> vector_dif
		// tento prvek/prvky pripoj k vektoru cells (pokud tam jiz neni) --> merge 

		vector_dif(neigh, neighi);		// ve vektorech neigh a neighi odebrany prvky, ktere se vyskytuji v obou vektorech
		merge(neigh, neighi);

		merge(cells, neigh);			// 
	} // END for (i; cells)


//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << cells[i] << " "; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 2) SECONDARY NEIGHBOURS
	// umet rozsirit list ovlivnenych bunek na list ovlivnenych paru bunek (tj. rozsirit o sekundarni sousedy), ...
	std::vector<int> sec;
	std::vector<int> sec_pos;
	std::vector<int> sio;
	std::vector<double> sap;
	std::vector<int> dif;

	count = 0; // vynulovani citace (novy vektor prislusny vektoru sec budeme cislovat opet od jednicky)
	//sec_neigh(cells, sec);
	//int kk = 0;
	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			//if (type == 2) { }
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede jsou vsichni v cells
		else {
			cell = con.compute_cell(c, cells_pos[2*i], cells_pos[2*i+1]);
			if (cell == true) {						// ______________________ sousede bunky v con, pokud existuji !!!
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {					// ______________________ sousede bunky v newcon, pokud existuji !!!
				c.neighbors(neighi);

			}

		} // END if..else (id)
		//kk++;
		
		dif.clear();
		vector_dif(neigh, cells, dif);
		merge(sec, dif);
		dif.clear();
		vector_dif(neighi, cells, dif);
		merge(sec, dif);
	} // END for (i; cells)

	// ulozeni pozic secondary particles
	for (i = 0; i < sec.size(); i++) {
		find_pos(ijk, q, sec[i], &con);		// pozice generatoru je stejna v con i v newcon
		sec_pos.push_back(ijk); sec_pos.push_back(q);
	}

	// !!!!!!! dulezity predpoklad: prvky vektoru sec se neshodji s prvky cells (tj. modifikovanymi casticemi)

//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
//		std::cout << sec[i]  << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// C) ENERGY

	// specifikovat energii a efektivne ji spocist pro bunky, pary bunek, ..., uvedene na poskytnutem liste
	// vypocet energie	- jednotlive bunky - loop pres cells
	//					- pary bunek - bunky v cells byli zmeneny, k nim potrebuji najit navic sekundarni sousedy, a pak loop pres pary; ale ne pres vsechny
	//							(neni potreba uvazovat pary sekundarni s.-sekundarni s.), bylo by vhodne tedy uvazovat dve struktury - cells a sekundarni s.
	//							a delat loop pres pary bunek z cells a takove, ze jedna bude z cells a druha ze sekundarnich
/*
// 1) FEASIBILITY
	if (feasibility(newcon, cells_pos, h_par[0], h_par[1], h_par[2])) {  }
	else { std::cout << "Unfeasible \n"; return 0; }
	if (h_par.size() > 3) {
		if (overlap_f(newcon, h_par[3])) {}
		else { std::cout << "Unfeasible (overlap) \n"; return 0; }
	}
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
// 2) CELL's ENERGY, energie pres jednotlive bunky
		// read histograms
	histogram hist;
	hist.read_hist();

	// sum up the energy
//double energy = 0;
	double val_bef, val_aft;
	double part1 = 0, part2 = 0, part3 = 0, part4 = 0, part5 = 0, part6 = 0; // ... // parts of energy
	//		volume		nof					 vol difs  				radius
	// energy signs > before the change = "-" , after the change = "+" 
	int j = 0;
	double val_a = 0, val_b = 0; // promenne pro souhrnne statistiky
	//std::vector<double> vec_a, vec_b;
	int tp = 0;
	// stejny cyklus jako pro zjisteni secondary neighbours
//# pragma omp parallel for shared(con, newcon)
	for (i = 0; i < cells.size(); i++) {
		val_bef = 0; val_aft = 0;
		if (cells[i] == id) {
			if (type == 1) {	// add
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					// 1. volume
					val_aft = c.volume();
					//part1 = part1 - val_aft;
					// 2. number of faces
					//val_aft = c.number_of_faces();
					//part2 = part2 - val_aft;
					// 6. radius
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

					val_a = val_a + val_aft;
					//vec_a.push_back(val_aft);
					tp++;
					info.hist_aft.hist_act(val_aft, 1);
//					std::cout << "con: O " << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////////////////////////////////////////
					//energy = energy - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 - V_function(val_aft, hist);
				}
			
			}
			if (type==2){		// delete
				j--;								// ______________________________________________________ inc or dec ????????????????
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					val_bef = c.volume();
					//part1 = part1 + val_bef;
					//val_bef = c.number_of_faces();
					//part2 = part2 + val_bef;
					//val_bef = con.p[ijk_mb][4 * q_mb + 3];

					val_b = val_b + val_bef;
					//vec_b.push_back(val_bef);
					tp--;
					info.hist_aft.hist_act(val_bef, 0);
//					std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: O \n"; ////////////////////////////////////////////////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 + V_function(val_bef, hist);
				}

			}
			if(type==3)	{		// move
				cell = con.compute_cell(c, ijk_mb, q_mb);
				if (cell == true) {
					val_bef = c.volume();
					//part1 = part1 + val_bef;
					//val_bef = c.number_of_faces();
					//part2 = part2 + val_bef;
					//val_bef = con.p[ijk_mb][4 * q_mb + 3];

					val_b = val_b + val_bef;
					//vec_b.push_back(val_bef);
					tp--;
					info.hist_aft.hist_act(val_bef, 0);
//					std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm); ////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1 + V_function(val_bef, hist);
					
				}
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					val_aft = c.volume();
					//part1 = part1 - val_aft;
					//val_aft = c.number_of_faces();
					//part2 = part2 - val_aft;
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

					val_a = val_a + val_aft;
					//vec_a.push_back(val_aft);
					tp++;
					info.hist_aft.hist_act(val_aft, 1);
//					std::cout << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) ; ////////////
					//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
					part1 = part1  - V_function(val_aft, hist);
				}
//				std::cout << "\n";

			}
		}	// end if (changed particle)					// zde nutne rozlisit pripady add/delete/move
		else {
			cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell == true) {	
				// 1. volume	
				val_bef = c.volume();
				//part1 = part1 + val_bef;
				// 2. number of faces
				//val_bef = c.number_of_faces();
				//part2 = part2 + val_bef;
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

				val_b = val_b + val_bef;
				//vec_b.push_back(val_bef);
				tp--;
				info.hist_aft.hist_act(val_bef, 0);

				part1 = part1 + V_function(val_bef, hist);
			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {	
				// 1. volume	
				val_aft = c.volume();
				//part1 = part1 - val_aft;
				// 2. number of faces
				//val_aft = c.number_of_faces();
				//part2 = part2 - val_aft;
				// 6. radius
				//val_aft = newcon.p[cells_pos[2 * j]][4 * cells_pos[2 * j + 1] + 3];

				val_a = val_a + val_aft;
				//vec_a.push_back(val_aft);
				tp++;
				info.hist_aft.hist_act(val_aft, 1);

				part1 = part1 - V_function(val_aft, hist);
			}
//			std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm)  << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " " << val_aft << " \n"; ////////////
			//energy = energy + sqrt(sqrt(sqrt(1 - hist.hist_value(val_bef) / hist.ocm))) - sqrt(sqrt(sqrt(1 - hist.hist_value(val_aft) / hist.ocm)));	// ____________________ znamenka
			//part1 = part1 + V_function(val_bef, hist) - V_function(val_aft, hist);
			// energy = energy + val_bef - val_aft;
			// energy = energy + min(val_bef, K) - min(val_aft, K);
			// energy = energy + (abs(val_bef - n0))/n0 - (abs(val_aft - n0))/n0;

			// energy = energy + V_function(val_bef) - V_function(val_aft);
		} // end else (changed particle)
		j++;
		//std::cout << "con: " << (1 - hist.hist_value(val_bef) / hist.ocm) << " & newcon: " << (1 - hist.hist_value(val_aft) / hist.ocm) << " \n"; ////////////
		//energy = energy + (1 - hist.hist_value(val_bef) / hist.ocm) - (1 - hist.hist_value(val_aft) / hist.ocm);	// ____________________ znamenka
	}

	
//	std::cout << "  energy (volume): " << part1 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	energy = energy + theta[1] * part1;
//	std::cout << "energy: " << energy << "\n";

//  COMPOUND ENERGY
	info.tp_aft = info.tp_bef + tp;
//	std::cout << info.tp_bef << " " << info.tp_aft << "\n";
	//std::cout << nonempty_cells(con) << " " << nonempty_cells(newcon) << "\n";
// -	energy = theta[0] * (V_function(info.tp_bef, hist) - V_function(info.tp_aft, hist));
// -	std::cout << energy << "\n";

	info.mean_aft[0] = info.mean_bef[0] - val_b + val_a;
//	std::cout << info.mean_bef[0] / info.tp_bef << " " << info.mean_aft[0] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[0] / nonempty_cells(con) << " " << info.mean_aft[0] / nonempty_cells(newcon) << "\n";
// -	energy = theta[1] * (V_function(info.mean_bef[0]/ info.tp_bef, hist) - V_function(info.mean_aft[0]/ info.tp_aft, hist));
//	std::cout << energy << "\n";


// -	info.varsum(newcon); // do info.var_aft ulozi novy rozptyl - nelze pocitat lokalne !!!  (...casove narocne...)
// -	energy = theta[1] * (V_function(info.var_bef[0]/ (info.tp_bef-1), hist) - V_function(info.var_aft[0]/ (info.tp_aft-1), hist));

	//std::cout << hist_dis(info.hist_bef, info.hist_aft) << "\n";
//	std::cout << hist_dis(info.hist_bef, hist) << " " << hist_dis(info.hist_aft, hist) << "\n";
// -	energy = theta[1] * (V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist));
//	std::cout << energy << "\n";

//	std::cout << "energy (nof): " << part2 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//energy = energy + theta[2] * part2;

	// ...

//	std::cout << "energy: " << energy << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// energie pro ruzne charakteristiky se uklada do ruznych promenych part1-6
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
// 3) IN/OUT INFORMATION - zastarale
	// when computing the energy over pairs, triplets, ... , one needs barycenter of union of cells to ensure the uniqueness of the contribution to the 
	//	energy; to compute the barycentrum the in/out information and the knowledge of the true coordinates (not periodic) is necessary
	//	therefore in the case of pair potential one needs in/out information for modified cells and secondary particles

	// algorithm:
	//		i)	zjistit informaci in/out ve stejnem okamziku kdy se vytvari vektory cells a sec (preferovana varianta)
	//		ii)	zjistit informaci in/out az pro dany seznam castic = cells+sec (viz nize)

	// entry: cells & sec vectors; the only particle we know the placement (IN) is the changed particle ID
//	int n = cells.size() + sec.size();	// number of particles whose placement has to be determined
//	n--;								// particle ID is IN

//	while (n > 0) {
//		n--;
//	}
*/
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// 4) PAIR POTENTIAL, energie pres pary, atd.
	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn;
		// dvojnasobny loop pres cells
	int k1 = 0, k2 = 0, ii = 0, jj = 0;
	int ijk2, q2;
	int j;
	bool cellc, celln, cell2c, cell2n;
	//int citac1 = 0, citac2 = 0;
	for (i = 0; i < cells.size(); i++) {
		k2 = 0; jj = 0;
		//k1 = i - citac1;
		//citac2 = 0;
		//ijk = cells_pos[2 * k1]; q = cells_pos[2 * k1 + 1];
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		if (cells[i] == id) {
			ii--;
			if (type == 1) {
				cellc = false;
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
			if (type == 2) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				celln = false;
				//citac1 = 1;
				//k1 = k1 - 1;

			}
			if (type == 3) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

		}

		for (j = 0; j < cells.size(); j++) {
			if (cells[i] < cells[j]) {			// prevents doublecounting   ............................................................ i & j
				//ijk2 = cells_pos[2 * k2]; q2 = cells_pos[2 * k2 + 1];
				//k2 = j - citac2;

				// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
				if (cells[j] == id) {
					if (type == 1) {
						cell2c = false;
						ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2, q2);

					}
					if (type == 2) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = false;
						//k2 = k2 - 1;
						//citac2 = 1;

					}
					if (type == 3) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = newcon.compute_cell(dn, cells_pos[2 * j], cells_pos[2 * j + 1]);

					}
				}
				else {
					ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);
					cell2n = newcon.compute_cell(dn, ijk2, q2);

				}

				// nyni mam dve bunky, pokud existuji, jak v con, tak i v newcon
				// staci overit, zda-li jsou sousede a pokud ano, tak spocist parovou energii

				if (cellc == true && cell2c == true) {				
					// v pripade delete/move neni pozice odstranovane castice ulozena v cells_pos, proto si ji musime ulozit zvlast - ijk2, q2
					if (are_neighbors(cc, ijk2, q2, &con)) {
						// compute pair energy
						//double V2(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz);
						energy = energy + V2_function(cc, dc);
						// !!! fce V2 vyzaduje SKUTECNE souradnice generatoru !!! ...........................................................................
					}
				}

				if (celln == true && cell2n == true) {
					if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
						// compute pair energy
						energy = energy - V2_function(cn, dn);   // jake znamenko???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);

			} // end..if(doublecounting)
			k2++; jj++;
		} // end..for(j; second loop over cells)
		k1++; ii++;
	} // end..for(i; first loop over cells)



		/*
				// pripad je-li jedna z nich id
				if (cells[k1] == id || cells[k2] == id) {
					// rozlis pripady add/delete/move
					if (type == 1) {	// add
						//cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]); // _________________ !!!!! ven do prvniho cyklu
						// takto to pocitam zbytecne mockrat, ale kdyz to dam ven musim pocitat tedy obe varianty (con, newcon) - cc, cn
						// a musim to ridit pres k1, k2
						//zajima me pouze con_new, a protoze pozice v cells_pos odpovidaji vektoru cells, pak nemusim rozlisovat, 
						// kt generator je id
						cell2 = newcon.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);

						if (celln == true && cell2 == true) {
							// compute pair energy - use only newcon
							//if are_neighbors
							//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
						}
					} // end..if (add)
					if (type == 2) {		// delete
						if (cells[i] == id) { 
							k1--;													// _________________________ k1--;
							cellc = con.compute_cell(cc, ijk_mb, q_mb);				// v cc neni spravna bunka
							cell2 = con.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
						}		
						if (cells[j] == id) { 
							k2--;													// _________________________ k2--;
							cell2 = con.compute_cell(d, ijk_mb, q_mb);				// v cc je spravna bunka
						}		
						
						if (cellc == true && cell2 == true) {
							// compute pair energy - use only con
						}
					} // end..if (delete)
					if (type == 3) {		// move
						//pozice v cells_pos neodpovidaji vektoru cells v pripade ze pracujeme s con
						//con
						if (cells[i] == id) {
							cellc = con.compute_cell(cc, ijk_mb, q_mb);
							cell2 = con.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);
							if (cellc == true && cell2 == true) {
								// compute pair energy - in con
							}
						}
						if (cells[j] == id) {
							cellc = con.compute_cell(cc, cells_pos[2 * i], cells_pos[2 * i + 1]);
							cell2 = con.compute_cell(d, ijk_mb, q_mb);
							if (cellc == true && cell2 == true) {
								// compute pair energy - in con
							}
						}
						//newcon
						cell2 = newcon.compute_cell(d, cells_pos[2 * j], cells_pos[2 * j + 1]);
						if (celln == true && cell2 == true) {
							// compute pair energy - in newcon
						}
					} // end..if (move)
					
				} // end if (cells[k1] == id || cells[k2] == id)
				// ani jedna neni id --> pozice bunek jsou zachovany, netreba rozlisovat operace
				else {
					//pouzivej k1, k2  (vsude)
					cell2 = con.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
					if (cellc == true && cell2 == true) {
						// compute pair energy - in con
					}
					cell2 = newcon.compute_cell(d, cells_pos[2 * k2], cells_pos[2 * k2 + 1]);
					if (celln == true && cell2 == true) {
						// compute pair energy - in newcon

						//if (are_neighbors(c, sr_add[2 * j], sr_add[2 * j + 1], &con)) {
						//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
					}

				} // end if..else (cells[k1] == id || cells[k2] == id)
				
			} // end..if (cells[k1] < cells[k2]) - against doublecounting
			k2++;
		} // end..for (second loop over cells)
		k1++;
	}
	*/
		
	// loop pres cells a sec
	// assumption: id does not neighbour with any particle from sec

	k1 = 0; ii = 0;
	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {
			ii--;
			if (type == 2) {
				k1--;
			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];

			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

			// mimo ID jsou souradnice ostatnich generatoru nezmeneny
//			if (io[ii] == 0) {
// //////	#pragma omp simd {    // - vektorizace - na co nejjednodussi operace                           ukazka VEKTORIZACE
//				x = con.p[ijk][4 * q];
//				y = con.p[ijk][4 * q + 1];
//				z = con.p[ijk][4 * q + 2];
// //////		}
//			}
//			else {
//				x = ap[3 * (io[ii] - 1)];
//				y = ap[3 * (io[ii] - 1) + 1];
//				z = ap[3 * (io[ii] - 1) + 2];
//			}
	
			for (j = 0; j < sec.size(); j++) {
				if (cells[i] < sec[j]) {			// prevents doublecounting
					ijk2 = sec_pos[2 * j]; q2 = sec_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);


					if (cellc == true && cell2c == true) {
						if (are_neighbors(cc, ijk2, q2, &con)) {
							// compute pair energy
							energy = energy + V2_function(cc, dc);
						}
					}
					cell2n = newcon.compute_cell(dn, ijk2, q2);

					if (celln == true && cell2n == true) {
						if (are_neighbors(cn, ijk2, q2, &newcon)) {
							// compute pair energy
							energy = energy - V2_function(cn, dn);
						}
					}
				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
		k1++; ii++;
	} // end..for (loop over cells)

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// -	energy = theta[0] * part1 + theta[1] * part2 + theta[2] * part3 + theta[3] * part4 + theta[4] * part5 + theta[5] * part6;
	// celkova energie je vysledkem slozeni nekolika casti doplnenych o vahy

	energy = theta[1] * energy;

	return exp(energy);
}




void LAG_container(container_poly &con, container_poly &newcon, int type, int id)
{
	//	[in,out]	con		container before the change
	//	[in,out]	newcon	container after the change (local variable, 
	//	[in]		type	suggested change (add/delete/move)
	//  [in]		id		change proposal (id of changed particle)
	
	// con je puvodni container

	// newcon je jeho kopie, v te provedeme zmenu (add/delete/move) - rekneme, ze zmena je jiz provedena
	// type reprezentuje zmenu: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)


	// 0) STRUCTURE OF CONS, structure must be kept identical (for both con and con_copy)
	// ADD - pridana castice je na konci nejakeho boxu v newcon, pozice ostatnich se shoduji
	// DELETE - castice se odstrani, pozice dalsich castic v boxu v newcon se posunou -> reseni: trik s odebranim a opetovnym pridanim
	if (type == 2) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}
	// MOVE - zmenena castice je v con i v newcon opet na ruzne pozici (move=delete+add), 
	// id castice bylo zachovano, ale pozice se zmeni, tim dojde k posunu pozic i ostatnich castic -> reseni: trik odeber a pridej
	if (type == 3) {
		int ijk_del, q_del;
		double x_del, y_del, z_del, rad_del;
		find_pos(ijk_del, q_del, id, &con);
		x_del = con.p[ijk_del][4 * q_del];      // coordinates of deleted particle
		y_del = con.p[ijk_del][4 * q_del + 1];
		z_del = con.p[ijk_del][4 * q_del + 2];
		rad_del = con.p[ijk_del][4 * q_del + 3];
		erase(ijk_del, q_del, &con);
		con.put(id, x_del, y_del, z_del, rad_del);
	}

}




bool LAG_cells(container_poly &con, container_poly &newcon, int type, int id, con_info &info, std::vector<int> &cells, std::vector<int> &cells_pos)
{
	//	[in]	con		container before the change
	//	[in]	newcon	container after the change (local variable, 
	//	[in]	type	suggested change (add/delete/move)
	//  [in]	id		change proposal (id of changed particle)
	//	[in]	info	structure containing summary characteristics of container
	//	[out]	cells
	//	[out]	cells_pos
	
	// con je puvodni container

	// newcon je jeho kopie, v te provedeme zmenu (add/delete/move) - rekneme, ze zmena je jiz provedena
	// type reprezentuje zmenu: 1-add, 2-delete, 3-move

	// porovnanim con a newcon urcime ovlivnene bunky - ziskame jejich id (z listu sousedu)
	// pro tyto bunky potrebujeme zjistit jejich ijk,q - nejlepe aby se shodovali  v obou kontejnerech ! - k tomu potrebujeme, aby menena castice
	//		byla na konci oddilu pole spravneho boxu (lze docilit trikem: smazat a znovu pridat) (pak budou mit ostatni castice stejnou polohu)


	// 1) LIST OF MODIFIED PARTICLES
	// list ovlivnenych bunek se muze lisit pro oba kontejnery (pripady add a delete) - poznacit si odlisnost? - odlisnost v jedne castici (te pridavane/odebirane)
	// zaklad listu lisicich se bunek tvori sousede zmenene castice

	int ijk, q;
	int ijk_ma, q_ma, ijk_mb, q_mb; // a - after (newcon), b - before (con)
	voronoicell_neighbor c;
	bool cell;

	//std::vector<int> empty; // empty cells
	info.empty.clear();


	// 1a) INITIAL PARTICLES
	if (type == 1) {
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell == true) {						//  muze se stat, ze pridany generator vytvori prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);

		} // otherwise mozaika se nezmenila = energie se nezmenila = return exp(0)
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
			return false;
		}
	}

	if (type == 2) {
		find_pos(ijk_mb, q_mb, id, &con);
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze odebirany generator tvoril prazdnou bunku !!! - mozaika se nezmeni
			c.neighbors(cells);

			// assumption: no empty cells in the previous configuration
		} // otherwise mozaika se nezmenila = energie se nezmenila = return 0
	}

	std::vector<int> neigh;
	bool cell2;

	if (type == 3) {
		find_pos(ijk_mb, q_mb, id, &con);				// pozice generatoru NEMUSI byt stejna v con i v newcon !!! (zmena souradnic muze zpusobit zmenu boxu)
		cell = con.compute_cell(c, ijk_mb, q_mb);
		if (cell == true) {						// muze se stat, ze generator tvoril prazdnou bunku !!! (assumption: there were no empty cells before)
			c.neighbors(cells);

		}
		find_pos(ijk_ma, q_ma, id, &newcon);
		cell2 = newcon.compute_cell(c, ijk_ma, q_ma);
		if (cell2 == true) {					//  muze se stat, ze zmeneny generator vytvori prazdnou bunku !!!
			c.neighbors(neigh);

		}
		else {
			info.empty.push_back(ijk_ma);
			info.empty.push_back(q_ma);
		}

		merge(cells, neigh);

	}


	//std::cout << "jsem tu \n";
	//	std::cout << cell << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	if (type == 3) { std::cout << cell2; } ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	for (int ii = 0; ii < cells.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << cells[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//	for (int ii = 0; ii < io.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << io[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	for (int ii = 0; ii < ap.size(); ii++) { ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << ap[ii] << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	
	// 1b) ALGORITHM for SEARCHING POSSIBLY MODIFIED PARTICLES
	// implementace algoritmu vyhledani zmenenych bunek
	// cells obsahuje ids bunek, ktere se nachazi jak v con tak i v newcon
	// funkce ktera porovna dva listy integeru a vrati jejich odlisnosti
	int i;
	std::vector<int> neighi, verti;
	
	bool shoda = false;
	// kontrolni mechanismus:
	int a, b;

	cells_pos.clear();

	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			if (type == 2) { cells_pos.push_back(-1); cells_pos.push_back(-1); } // delete (castice neni v newcon)
			else { cells_pos.push_back(ijk_ma); cells_pos.push_back(q_ma); } // uloz polohu z newcon
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede uz byli uvazovani
		else {
			find_pos(ijk, q, cells[i], &con);		// pozice generatoru je stejna v con i v newcon 
													//std::cout << "con: " << ijk << " " << q << " "; //////////////////////////////////////////////////////////////////////////////////////////////////
			cells_pos.push_back(ijk); cells_pos.push_back(q);
			cell = con.compute_cell(c, ijk, q);
			if (cell == true) {						// ______________________ sousede bunky id v con, pokud existuji !!! 
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, ijk, q);
			if (cell2 == true) {					// ______________________ sousede bunky id v newcon, pokud existuji !!! 
				c.neighbors(neighi);

			}
			else {
				info.empty.push_back(ijk);
				info.empty.push_back(q);
			}

			// kontrola:
			find_pos(a, b, cells[i], &newcon);		// pozice generatoru je stejna v con i v newcon ////////////////////////////////////////////////////////////
													//std::cout << "& newcon: " << a << " " << b << " \n"; /////////////////////////////////////////////////////////////////////////////////////////////
			if (a == ijk && b == q) {}
			else { std::cout << "ERROR: (LAG_recompute) not corresponding placement! \n"; } //////////////////////////////////////////////////////////


		} // END if..else (id)
		  // najdi v cem se lisi neigh a neighi --> vector_dif
		  // tento prvek/prvky pripoj k vektoru cells (pokud tam jiz neni) --> merge 

		vector_dif(neigh, neighi);		// ve vektorech neigh a neighi odebrany prvky, ktere se vyskytuji v obou vektorech
		merge(neigh, neighi);

		merge(cells, neigh);			// 
	} // END for (i; cells)


	  //	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //		std::cout << cells[i] << " "; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return true;
}



void LAG_sec(container_poly &con, container_poly &newcon, int id, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> &sec, std::vector<int> &sec_pos) 
{

	//	[in]	con		container before the change
	//	[in]	newcon	container after the change (local variable, 
	//  [in]	id		change proposal (id of changed particle)
	//	[in]	cells	modified cells
	//	[in]	cells_pos
	//	[out]	sec		secondary particles
	//	[out]	sec_pos


	// 2) SECONDARY NEIGHBOURS
	// umet rozsirit list ovlivnenych bunek na list ovlivnenych paru bunek (tj. rozsirit o sekundarni sousedy), ...
	//std::vector<int> sec;
	//std::vector<int> sec_pos;
	std::vector<int> dif;
	std::vector<int> neigh, neighi;

	int i, ijk, q;
	bool cell, cell2;

	voronoicell_neighbor c;

	
	for (i = 0; i < cells.size(); i++) {
		neigh.clear();
		neighi.clear();
		if (cells[i] == id) {
			//if (type == 2) { }
		}						// tato castice nemusi byt pritomna v obou kontejnerech, navic jeji sousede jsou vsichni v cells
		else {
			cell = con.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell == true) {						// ______________________ sousede bunky v con, pokud existuji !!!
				c.neighbors(neigh);

			}
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
			if (cell2 == true) {					// ______________________ sousede bunky v newcon, pokud existuji !!!
				c.neighbors(neighi);

			}

		} // END if..else (id)
		  //kk++;

		dif.clear();
		vector_dif(neigh, cells, dif);
		merge(sec, dif);
		dif.clear();
		vector_dif(neighi, cells, dif);
		merge(sec, dif);
	} // END for (i; cells)

	// ulozeni pozic secondary particles
	for (i = 0; i < sec.size(); i++) {
		find_pos(ijk, q, sec[i], &con);		// pozice generatoru je stejna v con i v newcon
		sec_pos.push_back(ijk); sec_pos.push_back(q);
	}

	// !!!!!!! dulezity predpoklad: prvky vektoru sec se neshodji s prvky cells (tj. modifikovanymi casticemi)

	//	for (int i = 0; i < cells.size(); i++) { //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		std::cout << sec[i]  << " "; ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	} ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	std::cout << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

}



bool LAG_feasibility(container_poly &con, container_poly &newcon, std::vector<double> &h_par, std::vector<int> cells_pos)
{

	//	[in]	con			container before the change
	//	[in]	newcon		container after the change (local variable, 
	//  [in]	id			change proposal (id of changed particle)
	//	[in]	h_par		hardcore parameters
	//	[in]	cells_pos	positions of modified particles



	// 3) ENERGY

	// specifikovat energii a efektivne ji spocist pro bunky, pary bunek, ..., uvedene na poskytnutem liste
	// vypocet energie	- jednotlive bunky - loop pres cells
	//					- pary bunek - bunky v cells byli zmeneny, k nim potrebuji najit navic sekundarni sousedy, a pak loop pres pary; ale ne pres vsechny
	//							(neni potreba uvazovat pary sekundarni s.-sekundarni s.), bylo by vhodne tedy uvazovat dve struktury - cells a sekundarni s.
	//							a delat loop pres pary bunek z cells a takove, ze jedna bude z cells a druha ze sekundarnich
	
	// 3a) FEASIBILITY
	if (feasibility(newcon, cells_pos, h_par[0], h_par[1], h_par[2])) {  }
//	else { std::cout << "Unfeasible \n"; return false; }
	else { return false; }
	if (h_par.size() > 3) {
		if (overlap_f(newcon, h_par[3])) {}
//		else { std::cout << "Unfeasible (overlap) \n"; return false; }
		else { return false; }
	}

	return true;
	
}





void LAG_V1(container_poly &con, container_poly &newcon, int type, int id, con_info &info, histogram &hist, histogram &hist2, std::vector<double> &parts, std::vector<int> cells, std::vector<int> cells_pos)
{

	//	[in]	con				container before the change
	//	[in]	newcon			container after the change (local variable, 
	//	[in]	type
	//	[in]	ijk_mb, q_mb 
	//  [in]	id				change proposal (id of changed particle)
	//	[in]	info
	//	[in]	hist, hist2
	//	[out]	parts			computed potentials
	//	[in]	cells		
	//	[in]	cells_pos		positions of modified particles

	// muzeme uvazovat vice potencialu (part1,...), souhrnnych statistik (val1_a, val1_b, ...), histogramu (hist1,...), V_funkci (V_function1,...)

	// 3b) CELL's ENERGY, energie pres jednotlive bunky

	long double val_bef, val_aft;
	double part1 = 0, part2 = 0, part3 = 0, part4 = 0, part5 = 0, part6 = 0; // ... // parts of energy (volume, nof, vol difs, radius, ... )
	// jednotlive parts jsou potencialy z kterych se sklada celkova funkce energie
	// energy signs > before the change = "-" , after the change = "+" 
	int i, ijk_mb, q_mb;
	bool cell, cell2;
	voronoicell_neighbor c;
	long double val1_a = 0, val1_b = 0, val2_a = 0, val2_b = 0; // promenne pro souhrnne statistiky
	//std::vector<double> vec_a, vec_b;
	int tp = 0;

	if (type == 1) {} else { find_pos(ijk_mb, q_mb, id, &con); }

	// stejny cyklus jako pro zjisteni secondary neighbours
	//# pragma omp parallel for shared(con, newcon)
	for (i = 0; i < cells.size(); i++) {
		val_bef = 0; val_aft = 0;
		if (cells[i] == id) {
			if (type == 1) {	// add
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
				if (cell == true) {
					tp++;										// aktualizace celkoveho poctu castic
					// 1. volume
					val_aft = c.volume();						// vypocet charakteristiky
					//val2_a = val2_a + val_aft;					// kumulovana suma pro aktualizaci vyberoveho prumeru
					info.hist2_aft.hist_act(val_aft, 1);			// aktualizace histogramu
//					part2 = part2 - V_function_2(val_aft, hist2);	// zmena potencialu
					// 2. number of faces
					val_aft = c.number_of_faces();
					val1_a = val1_a + val_aft;
					info.hist_aft.hist_act(val_aft, 1);
//					part1 = part1 - V_function(val_aft, hist);
					// ...

					// 6. radius
					//val_aft = newcon.p[cells_pos[2 * i]][4 * cells_pos[2 * i + 1] + 3];

				}

			}
			if (type == 2) {		// delete		
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

					// ...
				}

			}
			if (type == 3) {		// move
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
				cell = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
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
				
			}
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
			cell2 = newcon.compute_cell(c, cells_pos[2 * i], cells_pos[2 * i + 1]);
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


	//	std::cout << "  energy (volume): " << part1 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// - energy = energy + theta[1] * part1;
	//	std::cout << "energy: " << energy << "\n";

	parts[0] = part1; parts[1] = part2; // ...

	//  3c) COMPOUND ENERGY

	//	I) total particles
	info.tp_aft = info.tp_bef + tp;
	//	std::cout << info.tp_bef << " " << info.tp_aft << "\n";
	//std::cout << nonempty_cells(con) << " " << nonempty_cells(newcon) << "\n";
	// -	energy = theta[0] * (V_function(static_cast<double>(info.tp_bef), hist) - V_function(static_cast<double>(info.tp_aft), hist));
	// -	parts[0] = (V_function(static_cast<double>(info.tp_bef), hist) - V_function(static_cast<double>(info.tp_aft), hist));
	// -	std::cout << energy << "\n";
	//std::cout << tp << " ; " << val1_b << " " << val1_a << " ; " << val2_b << " " << val2_a << " \n";

	//	II) mean value
	info.mean_aft[0] = info.mean_bef[0] - val1_b + val1_a;
	//info.mean_aft[1] = info.mean_bef[1] - val2_b + val2_a;
	//info.mean_aft[1] = 1;
	//std::cout << info.mean_bef[0] / info.tp_bef << " " << info.mean_aft[0] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[1] / info.tp_bef << " " << info.mean_aft[1] / info.tp_aft << "\n";
	//std::cout << info.mean_bef[0] / nonempty_cells(con) << " " << info.mean_aft[0] / nonempty_cells(newcon) << "\n"; 
// -	parts[0] = (V_function(info.mean_bef[0]/ static_cast<double>(info.tp_bef), hist) - V_function(info.mean_aft[0]/ static_cast<double>(info.tp_aft), hist));		// nof
// -	parts[1] = (V_function_2(info.mean_bef[1]/ static_cast<double>(info.tp_bef), hist2) - V_function_2(info.mean_aft[1]/ static_cast<double>(info.tp_aft), hist2)); // vol
	double okno = 1;
// -	parts[1] = (V_function_2(okno / static_cast<double>(info.tp_bef), hist2) - V_function_2(okno / static_cast<double>(info.tp_aft), hist2)); // vol
	//	std::cout << energy << "\n";
	
	//	III) variance
	// -	info.varsum(newcon); // do info.var_aft ulozi novy rozptyl - nelze pocitat lokalne !!!  (...casove narocne...)
	// -	parts[0] = (V_function(info.var_bef[0]/ (info.tp_bef-1), hist) - V_function(info.var_aft[0]/ (info.tp_aft-1), hist));

	//	IV) histogram discrepancy
	//std::cout << hist_dis(info.hist_bef, info.hist_aft) << "\n";
	//	std::cout << hist_dis(info.hist_bef, hist) << " " << hist_dis(info.hist_aft, hist) << "\n";
	parts[0] = (V_function(info.hist_bef, hist) - V_function(info.hist_aft, hist));
	parts[1] = (V_function_2(info.hist2_bef, hist2) - V_function_2(info.hist2_aft, hist2));
	//	std::cout << energy << "\n";

	//	std::cout << "energy (nof): " << part2 << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//energy = energy + theta[2] * part2;

	// ...

	//	std::cout << "energy: " << energy << "\n"; ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// energie pro ruzne charakteristiky se uklada do ruznych promenych part1-6

}




double LAG_V2(container_poly &con, container_poly &newcon, int type, int id, con_info &info, std::vector<int> cells, std::vector<int> cells_pos, std::vector<int> sec, std::vector<int> sec_pos)
{

	//	[in]	con				container before the change
	//	[in]	newcon			container after the change (local variable, 
	//	[in]	type
	//	[in]	ijk_mb, q_mb
	//  [in]	id				change proposal (id of changed particle)
	//	[in]	h_par			hardcore parameters
	//	[in]	cells_pos		positions of modified particles

	// 3c) PAIR POTENTIAL, energie pres pary, atd.
	double energy = 0;

	voronoicell_neighbor dc, dn, cc, cn;
	double xc, yc, zc, xn, yn, zn, xxc, yyc, zzc, xxn, yyn, zzn;
	// dvojnasobny loop pres cells
	int ijk, q, ijk2, q2, ijk_mb, q_mb;
	int i, j;
	bool cellc, celln, cell2c, cell2n, arnc, arnn;
	int citac1 = 0, citac2 = 0;

	if (type == 1) {} else { find_pos(ijk_mb, q_mb, id, &con); }

//	std::cout << " cells+cells particles \n";
	for (i = 0; i < cells.size(); i++) {
		
		//ijk = cells_pos[2 * k1]; q = cells_pos[2 * k1 + 1];
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		if (cells[i] == id) {

			if (type == 1) {
				cellc = false;
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
			if (type == 2) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				celln = false;
			
			}
			if (type == 3) {
				cellc = con.compute_cell(cc, ijk_mb, q_mb);
				ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
				celln = newcon.compute_cell(cn, ijk, q);

			}
		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];
			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

		}

		for (j = 0; j < cells.size(); j++) {
			if (cells[i] < cells[j]) {			// prevents doublecounting   ............................................................ i & j
												//ijk2 = cells_pos[2 * k2]; q2 = cells_pos[2 * k2 + 1];
				
				// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
				if (cells[j] == id) {
					if (type == 1) {
						cell2c = false;
						ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
						cell2n = newcon.compute_cell(dn, ijk2, q2);

					}
					if (type == 2) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = false;
					
					}
					if (type == 3) {
						cell2c = con.compute_cell(dc, ijk_mb, q_mb);
						ijk2 = ijk_mb; q2 = q_mb;
						cell2n = newcon.compute_cell(dn, cells_pos[2 * j], cells_pos[2 * j + 1]);

					}
				}
				else {
					ijk2 = cells_pos[2 * j]; q2 = cells_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);
					cell2n = newcon.compute_cell(dn, ijk2, q2);

				}

				// nyni mam dve bunky, pokud existuji, jak v con, tak i v newcon
				// staci overit, zda-li jsou sousede a pokud ano, tak spocist parovou energii

				arnc = false;
				if (cellc == true && cell2c == true) {
					// v pripade delete/move neni pozice odstranovane castice ulozena v cells_pos, proto si ji musime ulozit zvlast - ijk2, q2
					arnc = are_neighbors(cc, ijk2, q2, &con);
					if (arnc) {
						// compute pair energy
						//double V2(voro::voronoicell_neighbor &rc1, voro::voronoicell_neighbor &rc2, double &rx, double &ry, double &rz, double &rxx, double &ryy, double &rzz);
						energy = energy + V2_function(cc, dc);
						// !!! fce V2 vyzaduje SKUTECNE souradnice generatoru !!! ...........................................................................
					}
				}

				arnn = false;
				if (celln == true && cell2n == true) {
					//if (are_neighbors(cn, cells_pos[2 * k2], cells_pos[2 * k2 + 1], &newcon)) {
					arnn = are_neighbors(cn, ijk2, q2, &newcon);
					if (arnn) {
						// compute pair energy
						energy = energy - V2_function(cn, dn);   // jake znamenko???!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					}
				}
				//energy = energy + V2(c, d, x, y, z, xn, yn, zn);
//				if (arnc || arnn) { std::cout << "\n"; citac1++; }
			} // end..if(doublecounting)
		} // end..for(j; second loop over cells)
	} // end..for(i; first loop over cells)
//	std::cout << "number of pairs: " << citac1 << " \n";

	// loop pres cells a sec
	// assumption: id does not neighbour with any particle from sec
//	std::cout << " cells+sec particles \n";

	for (i = 0; i < cells.size(); i++) {
		// castice ID nemusi byt pritomna v obou con i newcon - add: chybi v con, delete: chybi v newcon
		// castice ID ale take nema sousedy mezi casticemi v sec = castici ID muzeme vynechat !!! 
		if (cells[i] == id) {

		}
		else {
			ijk = cells_pos[2 * i]; q = cells_pos[2 * i + 1];

			cellc = con.compute_cell(cc, ijk, q);
			celln = newcon.compute_cell(cn, ijk, q);

			// mimo ID jsou souradnice ostatnich generatoru nezmeneny
			//			if (io[ii] == 0) {
			// //////	#pragma omp simd {    // - vektorizace - na co nejjednodussi operace                           ukazka VEKTORIZACE
			//				x = con.p[ijk][4 * q];
			//				y = con.p[ijk][4 * q + 1];
			//				z = con.p[ijk][4 * q + 2];
			// //////		}
			//			}
			//			else {
			//				x = ap[3 * (io[ii] - 1)];
			//				y = ap[3 * (io[ii] - 1) + 1];
			//				z = ap[3 * (io[ii] - 1) + 2];
			//			}

			for (j = 0; j < sec.size(); j++) {
// 				if (cells[i] < sec[j]) {			// prevents doublecounting
					ijk2 = sec_pos[2 * j]; q2 = sec_pos[2 * j + 1];
					cell2c = con.compute_cell(dc, ijk2, q2);

					arnc = false;
					if (cellc == true && cell2c == true) {
						arnc = are_neighbors(cc, ijk2, q2, &con);
						if (arnc) {
							// compute pair energy
							energy = energy + V2_function(cc, dc);
						}
					}
					cell2n = newcon.compute_cell(dn, ijk2, q2);

					arnn = false;
					if (celln == true && cell2n == true) {
						arnn = are_neighbors(cn, ijk2, q2, &newcon);
						if (arnn) {
							// compute pair energy
							energy = energy - V2_function(cn, dn);
						}
					}
//					if (arnc || arnn) { std::cout << "\n"; citac2++; }
//				} // end..if(doublecounting)
			} // end..for(second loop = loop over sec)
		} // end..if..else (cells[i] = id)
	} // end..for (loop over cells)
//	std::cout << "number of pairs: " << citac2 << " \n";

	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	  // -	energy = theta[0] * part1 + theta[1] * part2 + theta[2] * part3 + theta[3] * part4 + theta[4] * part5 + theta[5] * part6;
	  // celkova energie je vysledkem slozeni nekolika casti doplnenych o vahy

	//energy = theta[1] * energy;

	//std::cout << energy << "\n";
	return energy; // returns potential
	
}





void LAG_estim(double &th_estim, double &z_estim, container_poly &con, container_poly &newcon, int N, std::vector<double> &hpar_estim, con_info &info, histogram &hist, histogram &hist2, double lx, double ux, double ly, double uy, double lz, double uz)
{
	// [out]	th_estim					the estimate of parameter theta
	// [out]	z_estim						the estimate of parameter zet
	// [in]		con							the container with stored particles containing the particles radii
	// [in]		newcon						copy of the container con
	// [in]		N							the size of the sample for estimating integrals by Monte Carlo
	// [in]		hpar_estim					the estimates of hardcore parameters
	// [in]		lb, ub						lower and upper bound of subwindow (different bounds in every coordinate)

	int id = 0;
	int r, i;
	int na = 0;
	double eloc_e;
	bool add;
	double x, y, z, rr;
	std::vector<int> va;
	std::vector<double> vb, vc;

	va.clear();
	vb.clear();
	vc.clear();

	// 0) ESTIM SUBWINDOW
	//estimates are done from a subwindow, i.e. removable points are computed on subwindow only
	// - bounds of subwindow are specified by doubles lb, ub (lb=0,ub=1 means to use the whole window)
	// - restriction on the subwindow helps to avoid boundary effects

	// 1A) no_removable returns number of removable points in container
	// - !! na rozdil od pridani bodu, ktere muze porusit horni mez, odebrani bodu dolni mez porusit nemuze
	//std::cout << "Number of removable: " << no_removable(rcon, 0, rbeta_e + 0.00001, rB_e, vb, 0, 1) << " from " << rcon.total_particles() << "\n";
	// - the question is whether use the correction rbeta_e + 0.00001 or not
	// - fc no_removable computes coefficients b(x)
//	r = no_removable(con, newcon, hpar_estim, vb, lx, ux, ly, uy, lz, uz);
//	r = nonempty_cells(con);
	r = LAG_removable(con, newcon, vb, hpar_estim, info, hist, hist2, lx, ux, ly, uy, lz, uz);
	// POTREBUJI 3 VEKTORY VB, nebo ukladat trojice do jednoho vektoru
	if (r == 0) { std::cout << "ERROR: Estimates can not be computed (r)! \n"; return; }
	else { std::cout << "Number of removable (subwindow): " << r << "\n"; }

	// 1B) vypocet integralu - koeficientu a(x), c(x) (PLL_coeff)
	for (i = 0; i < N; i++) {
		x = lx + (ux - lx)*uniform(0, 1); y = ly + (uy - ly)*uniform(0, 1); z = lz + (uz - lz)*uniform(0, 1);
		//rr = triangle(0.005, 0.03, 0.02); // generating of the sample - coordinates from (lb,ub) and radius from ... 
		//rr = 0.005 + 0.025*uniform(0, 1); 
		rr = 0.05*uniform(0, 1);
		//rr = gamma(3, 0.5);
		// or using the average radius in container:
		//rr = ave_rad(conp);
		//rr = 100*uniform(0,1); 
		// (the same sampler for radius as in try_add fc is used)

//		add = addable(rcon, x, y, z, rr, ralfa_e, rbeta_e, rB_e, riota_e); // try if the particle is addable or not 
		//LAG_container(con, newcon, 1, id);
		newcon.put(id, x, y, z, rr);
		std::vector<int> cells; std::vector<int> cells_pos;
		LAG_cells(con, newcon, 1, id, info, cells, cells_pos);

		std::vector<int> sec; std::vector<int> sec_pos;
		LAG_sec(con, newcon, id, cells, cells_pos, sec, sec_pos);
		bool pripustnost = true;
		//pripustnost = LAG_feasibility(conp, conp_copy, hard_par, cells_pos); 
		if (pripustnost == false) { 
			// neni addable
			va.push_back(0); // save this information
			// loc energy is infite = point is not addable
			vc.push_back(0); // save arbitrary value - in this case the value is not important
			vc.push_back(0);
			vc.push_back(0);
		}
		else {
			va.push_back(1); // save this information
			std::vector<double> parts; parts.resize(3); // parts je nynu vektor tri nul
			LAG_V1(con, newcon, 1, id, info, hist, hist2, parts, cells, cells_pos); // hist se p
			//pst = exp(theta[0] * parts[0] + theta[1] * parts[1]);
			//vc.push_back(-log(eloc_e)); // save the correct value
			parts[2] = LAG_V2(con, newcon, 1, id, info, cells, cells_pos, sec, sec_pos);
			// POTREBUJI 3 VEKTORY VC, nebo muzu ukladat trojice do jednoho vektoru
			vc.push_back(parts[0]); // save local energies
			vc.push_back(parts[1]);
			vc.push_back(parts[2]);
		}
		//add = 0;
		//va.push_back(add); // save this information
		//
		//if (add == 0) {
		//	// loc energy is infite = point is not addable
		//	vc.push_back(0); // save arbitrary value - in this case the value is not important
		//}
		//else {
//		//	eloc_e = try_add(0, x, y, z, rr, rcon, 1, ralfa_e, rbeta_e, rB_e, riota_e); // compute the adding energy iff addable
		//	eloc_e = 1;
		//	vc.push_back(-log(eloc_e)); // save the correct value
		//}
	}

	for (i = 0; i < N; i++) { na = na + va[i]; }
	if (na == 0) { std::cout << "ERROR: Estimates can not be computed (a)! \n"; return; }
	std::cout << "Number of addable: " << na << " from " << N << "\n";

	std::cout << "\n";

	// 2) Takacz-Fiksel estimate can be obtained for different choices of testing functions h1 and h2
	// 2A)------------------------------------------------------ Quick estimate --------------------------------------------------------------------------
	// - first testing fc is identically equal to 1
	//bisection(va, vb, vc, k, 0);
	//secant(va, vb, vc, k, 0);
/*	th_estim = NR(va, vb, vc, r, N);		// 3) Newton-Raphson method
	 
	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(-th_estim*vc[i])*va[i];
	}
	z_estim = (N*r) / z;
	// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

	std::cout << "Q: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
*/

	// 2B)---------------------------------------------------- Max. pseudolikelihood estimate ------------------------------------------------------------
	// 2B.I) classic approach searching for roots of derivation equation
//	std::vector<double> wc;
//	int siz = vc.size() / 3;
//	wc.resize(siz);
//	for (i = 0; i < siz; i++) { wc[i] = vc[3 * i + 2]; }

	// assumption: zajima nas pouze odhad theta2, tj koeficientu u part2, ostatni thety uvazujme jako zname konstanty
	double th1_const = 1, th2_const = 1;

	//bisection(va, vb, vc, k, 0);							// 1) bisection method
	//secant(va, vb, vc, k, 0);								// 2) secant method
	th_estim = NR(va, vb, vc, r, 0);						// 3) Newton-Raphson method

	z = 0;
	for (i = 0; i < N; i++) {
		z = z + exp(th1_const*vc[3*i] + th2_const*vc[3*i+1] + th_estim*vc[3*i+2])*va[i];
	}
	z_estim = (N*r) / z;
	// zet_estim(N, th_estim, va, vc); // estimation of intesnsity parameter --> testing fc is identically equal 1

	std::cout << "MPLE-c: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";

	// 2B.II) approach using LBFGS optimization algorithm directly to pseudolikelihood fc insted of its derivatives
	// - uses the library LBFGS++ described here: https://yixuan.cos.name/LBFGSpp/	
	/*
	//	0) LBFGS method
	const int n = 10;
	// Set up parameters
	LBFGSpp::LBFGSParam<double> param;
	param.epsilon = 1e-6;
	param.max_iterations = 100;

	// Create solver and function object
	LBFGSpp::LBFGSSolver<double> solver(param);
	PLL fun(n, va, vb, vc, r);

	// Initial guess
	Eigen::VectorXd ex = Eigen::VectorXd::Zero(n);
	// x will be overwritten to be the best point found
	double fx;
	int niter = solver.minimize(fun, ex, fx);

	std::cout << niter << " iterations" << std::endl;
	std::cout << " x = \n" << ex.transpose() << std::endl;
	std::cout << " f(x) = " << fx << std::endl;

	// ?? x contains estimates of theta and zet --> how?
	//th_estim = x[0]; z_estim = x[1];
	//std::cout << "MPLE-LBFGS: Estimates of (theta,z): " << th_estim << " " << z_estim << "\n";
	//std::cout << "            Number of iterations was " << niter << " and optimal value reached was " << fx << "\n";
	*/


	// 3) using an initial guess on parameters (Takacz-Fiksel estimate can be used) the maximum likelihood procedure can be started
	// - improves the estimate values from the initial guess
	// - function MLE under construction below


}


int LAG_removable(container_poly &con, container_poly &con_copy, std::vector<double> &vb, std::vector<double> &hpar_estim, con_info &info, histogram &hist, histogram &hist2, double lx, double ux, double ly, double uy, double lz, double uz)
{
	//r = no_removable(con, newcon, hpar_estim, vb, lx, ux, ly, uy, lz, uz);

	int no = 0;
	int j, i, ci, id;
	int ijk, q;
	double x, y, z, r;

	for (j = 0; j < con.nxyz; j++) { // loop over boxes
		std::cout << j << " : " << con.co[j] << "\n";
		ci = 0; // poradove cislo castice; indikator poctu neodstranitelnych bodu = pouze pro kontrolu, zda-li jde o doplnek do celkoveho poctu
		for (i = 0; i < con.co[j]; i++) { // loop over particles in considered box - prtz se castice smaze a pak zase prida,
										   // je potreba co[j]-krat vzit prvni castici v bodu, ta se vzdy presune na konec
			id = con.id[j][ci];
			// std::cout << rcon.co[j] << " "; ///////////////////////////////////////////////////////////////////////////
			//std::cout << id << " ";  /////////////////////////////////////////////////////////////////////////////////
			//if (removable(rcon, id, ralfa, rbeta, rB)) { 
			find_pos(ijk, q, id, &con);	// find position of the particle
			// ijk = j, q = ci

			x = con.p[ijk][4 * q]; y = con.p[ijk][4 * q + 1]; z = con.p[ijk][4 * q + 2]; r = con.p[ijk][4 * q + 3]; // urci jeji souradnice a polomer
																														//		std::cout << "DELETE " << del << " " << conp.p[ijk][4 * q] << " " << conp.p[ijk][4 * q + 1] << " " << conp.p[ijk][4 * q + 2] << " " << conp.p[ijk][4 * q + 3] << "\n"; // popis cinosti

			//no = conp.total_particles();								// pocet castic v kontejneru \neq pocet bunek mozaiky !!!
																		// no = empty_cells(conp);

			erase(ijk, q, &con_copy);									// smazani v kopii, id ve fid neni potreba uvolnovat

			//		pst = LAG_recompute(conp, conp_copy, typ, del, hard_par, theta, info);	// prepocet energie
			//		pst = try_delete(ijk, q, conp, theta, alfa, beta, B, iota);	// spocte pravdepodobnost se kterou dojde k operaci DELETE
			// fci LAG_recompute si lze "poskladat":
			LAG_container(con, con_copy, 2, id);
			std::vector<int> cells; std::vector<int> cells_pos;
			LAG_cells(con, con_copy, 2, id, info, cells, cells_pos);
			std::vector<int> sec; std::vector<int> sec_pos;
			LAG_sec(con, con_copy, id, cells, cells_pos, sec, sec_pos);
			bool pripustnost = true;

			if (con.p[ijk][4 * q] > lx && con.p[ijk][4 * q] < ux && con.p[ijk][4 * q + 1] > ly && con.p[ijk][4 * q + 1] < uy && con.p[ijk][4 * q + 2] > lz && con.p[ijk][4 * q + 2] < uz) {
				// if the coordinates are in subwindow
			pripustnost = LAG_feasibility(con, con_copy, hpar_estim, cells_pos);
				if (pripustnost == false) {}
				else {
					std::vector<double> parts; parts.resize(3);
					LAG_V1(con, con_copy, 2, id, info, hist, hist2, parts, cells, cells_pos); // hist se p
					//pst = exp(theta[0] * parts[0]); // +theta[1] * parts[1]);
					parts[2] = LAG_V2(con, con_copy, 2, id, info, cells, cells_pos, sec, sec_pos);
					//pst = pst*exp(theta[2] * parts[2]);
					//std::cout << parts[0] << " " << parts[1] << " " << parts[2] << "\n";

//					vb.push_back(parts[0]); 
//					vb.push_back(parts[1]);
					vb.push_back(parts[2]); // save this information
					no++;
				}

			}
			else { ci++; }
			//ci++;
			//}
			con_copy.put(id, x, y, z, r);
		}
	}
	//std::cout << "\n"; ////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout << ci << "\n";
	return no;
}