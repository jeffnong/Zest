
#include <string>
#include <vector>
#include <map>

#pragma once

using namespace std;

namespace optics {

	typedef struct WaveLength {
		double lambda;
		string designation;
		string desc;
		string element;
		string source;
	} WL;

	static WL WLFQ[] = {
		{ 2.32542,"2.32542","infrared mercury line","Hg","" },
		{ 1.97009,"1.97009","infrared mercury line","Hg","" },
		{ 1.529582, "1.529582", "infrared mercury line","Hg","" },
		{ 1.060, "1.060" ,"Nd glass laser","Nd","Nd Laser" },
		{ 1.01398, "t" ,"infrared mercury line","Hg","" },
		{ 0.85211, "s" ,"infrared cesium line","Hd","" },
		{ 0.7065188, "r" ,"infrared cesium line","Cs","" },
		{ 0.6562725, "C" ,"red hydrogen line","H","" },
		{ 0.6438469, "C'" ,"red cadmium line","Cd","" },
		{ 0.6328, "HeNe","helium-neon-gas-laser","He-Ne","" },
		{ 0.5892938, "D","center of double sodium line","Na","" },
		{ 0.5875618, "d","yellow helium line","He","" },
		{ 0.546074, "e","green mercury line","Hg","" },
		{ 0.4861327, "F","blue hydrogen line","Hg","" },
		{ 0.4799914, "F'","blue cadmium line","Cd","" },
		{ 0.4358343, "g","blue mercury line","Hg","" },
		{ 0.4046561, "h","violet mercury line","Hg","" },
		{ 0.365, "i","utraviolet mercury line","Hg","" },
		{ 0.3341478, "0.3341478","utraviolet mercury line","Hg","" },
		{ 0.3125663, "0.3125663","utraviolet mercury line","Hg","" },
		{ 0.2967278, "0.2967278","utraviolet mercury line","Hg","" },
		{ 0.2804, "0.2804","ultraviolet mercury line","Hg","" },
		{ 0.2483, "0.2483","ultraviolet mercury line","Hg","" },
		{ 0.488 , "Argon" ,"","","" },
		{ 0.5145, "Argon 2" ,"","","" },
		{ 1.0641 , "ND:YAG" ,"","","" },
		{ 1.054 , "ND:Glass'" ,"","","" },
		{ 0.760 , "Ti:Al203" ,"","","" },
		{ 0.6943 , "Ruby" ,"","","" }
	};

	static std::map<string, WL> wlMap;

	class Glass
	{
	public:
		Glass(string Name, string Catalog); // Construct a glass
		double Index(double lambda); // Calculate Index from cofficients and fomula
		double Index(string line);
		double Pxy(string x, string y);
		double dPgF();

		string comment;
		string name;
		string catalog;
		int gid;
		string milnum;
		int ExeSub, status, meltfreq;

		double nd, nf, nc;
		double vd, xPgF, dxPgF;
		double density;
		double RelCost, CR, FR, SR, AR, PR;
		double lmin, lmax;
		enum FORMULA { schott = 1, sellmeier1 = 2, herzberger = 3, sellmeier2 = 4, Conrady = 5, sellmeier3 = 6, handbookoptics1 = 7, handbookoptics2 = 8, sellmeier4 = 9, extend = 10, sellmeier5 = 11, extend2 = 12 };
		enum StatusCode { standard = 0, perfer = 1, obsolete = 2, special = 3, melt = 4 };
		int IngExp, ExcSub;

		int Formula;

		double C[10];
		string CoffName[10];

		double D0, D1, D2, E0, E1, Ltk, TCE, TCE2, Temp;

		double T[21], D[21], L[21];

		double Schott(double lambda);
		double Sellmeier1(double lambda);
		double Herzberger(double lambda);
	};

	class GlassCat {

		friend class Glass;
		friend class Lens;
	public:
		static const string CatPath;
		GlassCat(string catname);
	public:
		string name;
		string comment;
		int number;
		vector<Glass> vecGlass;
		static void LoadGlassCat(GlassCat& cat);
		Glass FindGlass(string name);
		static Glass FindGlass(string cat,string glass);

		void ExportCSV(string path);
		void ExportJson(string path);		
	};
}