#include <stdio.h>
#include <fstream> // for ifstream
#include <time.h>
#include <chrono>
#include <random>
#include <iostream>
#include <string>
#include "Linux/Release/Tessellation.h"
using namespace std;

namespace GEN {
	uint64_t timeSeed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
	std::seed_seq ss{ uint32_t(timeSeed & 0xffffffff), uint32_t(timeSeed >> 32) };
	std::mt19937_64 rng(ss);
	std::uniform_real_distribution<double> unif(0, 1);
	std::uniform_real_distribution<double> stredy(0, 0.1);
} using namespace GEN;
// ^ set up the random number generator

namespace Params {
	string folder_out = ".";
		// output folder
	string name_out = "test";
		// name of output file, without ext
}

void Test(int begin_number, int iterations, bool remove = true, string out_path = "") {
	Tessellation Tessel;
	//Tessel.PrintEdges(SekceParams::RelPathData + L"/" + SekceParams::SouborVystup + to_wstring((*ind)++) + L".txt");

	for (int i = -1; ++i < begin_number; ) {
		Tessel.InsertGrain(unif(rng), unif(rng), unif(rng), 0.1);
	}
	clock_t begin = clock();
	cout << "Size " << Tessel.Grains.size() << ", press enter to start testing.";
	getchar();
	for (int i = -1; ++i < iterations; ) {
		Tessel.InsertGrain(unif(rng), unif(rng), unif(rng), 0.1);
		if(remove == true) Tessel.RemoveGrain(Tessel.Grains.back());
		if ((i + 1) % 10 == 0) cout << "step " << i + 1 << " size " << Tessel.Grains.size() << endl;
	};
	cout << "end, elapsed " << double(clock() - begin) / CLOCKS_PER_SEC << endl;
	if (out_path != "") {
		Tessel.PrintEdges(Params::folder_out + "/" + out_path + ".txt");
		cout << "saved " << Params::name_out << endl;
	}
}


class MojeTess : public Tessellation {

};

class MojeGrain : public Grain {
public: int DejPocetHran() {
	//...
}
		using Grain::Grain;
};

int main(int argc, char* argv[])
{
	Test(100, 100, true, "test"); // vysledek 24.165 sec, 1000 pridani, bez odebirani
	//return 0;
	// ^ use this to test the library

	Tessellation Tessel;
	int ret;
	for (int i = 0; i < 10; i++) {
		ret = Tessel.InsertGrain(GEN::unif(rng), GEN::unif(rng), GEN::unif(rng), GEN::stredy(rng));
		if (ret == 0) *Tessel.Grains.back() = *Tessel.Grains.back()->to_connected();
		// ^ We demonstrate here it is safe to add new grains when some grains in the tessellation have been transformed to te connected form.
	}
	Tessel.PrintEdges(Params::folder_out + "/" + Params::name_out+"10con.txt", true);
	// ^ print using connected interpretation ; setting only_once = true has no effect unless some grains have connected form
	Tessel.PrintEdges(Params::folder_out + "/" + Params::name_out + "10uncon.txt", false);
	// ^ print using unconnected interpretation

	Tessel.InsertGrain(shared_ptr<MojeGrain>(new MojeGrain(unif(rng), unif(rng), unif(rng), stredy(rng))));
	//  ^ this is how you insert customized grain

    return 0;
}
