
#include "stdafx.h"
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <vector>
#include <list>
#include "Common_Routine.h"
#include "glass.h"
#include "lens.h"

using namespace std;
using namespace Common;

namespace optics {

	const string GlassCat::CatPath = "C:\\Users\\Jeff\\Documents\\ZEST\\Glasscat\\";

	Glass::Glass(string _name, string _catalog) :name(_name), catalog(_catalog) {

	}

	double Glass::Index(double lambda)
	{
		double index;

		switch (Formula)
		{
		case 1:
			index = Schott(lambda);
			break;
		case 2:
			index = Sellmeier1(lambda);
			break;
		case 3:
			index = Herzberger(lambda);
			break;
		default:
			index = 1.0;
		}
		return index;

	}

	double Glass::Index(string line) {
		double l = wlMap[line].lambda;
		return Index(l);
	}

	double Glass::Pxy(string x, string y) {
		double nx = Index(x);
		double ny = Index(y);
		static double nf = Index("F");
		static double nc = Index("C");
		double Pxy = (nx - ny) / (nf - nc);
		return Pxy;
	}
	double Glass::dPgF() {
		//double dPCt = Pxy("C", "t") - (0.5450 + 0.004743*vd);
		//double dPCs = Pxy("C", "s") - (0.4029 + 0.002331*vd);
		//double dPFe = Pxy("F", "e") - (0.4884 - 0.000526*vd);
		//double dPig = Pxy("i", "g") - (1.7241 - 0.008382*vd);
		double dPgF = Pxy("g", "F") - (0.6438 - 0.001628*vd);

		return dPgF;
	}
	double Glass::Schott(double lambda)
	{
		double l2 = lambda*lambda;
		double l4 = l2*l2;
		double l6 = l2*l4;
		double l8 = l4*l4;
		double n2 = C[0] + C[1] * l2 + C[2] / l2 + C[3] / l4 + C[4] / l6 + C[5] / l8;
		return sqrt(n2);
	}
	double Glass::Sellmeier1(double lambda)
	{
		double l2 = lambda*lambda;
		double n2 = 1 + C[0] * l2 / (l2 - C[1]) + C[2] * l2 / (l2 - C[3]) + C[4] * l2 / (l2 - C[5]);

		return sqrt(n2);
	}

	double Glass::Herzberger(double lambda)
	{
		double L = 1 / (lambda*lambda - 0.028);
		double l2 = lambda*lambda;
		double n = C[0] + C[1] * L + C[2] * L*L + C[3] * l2 + C[4] * pow(lambda, 4) + C[5] * pow(lambda, 6);
		return n;
	}
	GlassCat::GlassCat(string _name) : name(_name) {
		if (vecGlass.size() == 0) {
			LoadGlassCat(*this);
		}
	}
	
	Glass GlassCat::FindGlass(string glass) {

		for (int i = 0; i < vecGlass.size(); i++) {
			if (vecGlass[i].name == glass) {
				Glass g = vecGlass[i];
				return g;
			}
		}
		return Glass("", "N.A.");
	}

	Glass GlassCat::FindGlass(string cat, string glass) {
		Glass g("", "N.A.");
		for (int i = 0; i < Lens::gcatList.size(); i++) {
			g = Lens::gcatList[i].FindGlass(glass);
			if (g.name != "")
				break;
		}
		return g;
	}

	void GlassCat::LoadGlassCat(GlassCat& cat)
	{
		string line;
		string mnemonic;

		string Name, catComment, gComment;
		int Formula;
		string MILID;
		double nd, vd;

		int excsub, status, meltfreq;
		double TCE, TCE2, density, dPgF;
		int IngExp;
		double C[10];
		double D0, D1, D2, E0, E1, Ltk, TempUsed;
		double RelCost, CR, FR, SR, AR, PR;
		double lmin, lmax;
		static const int N_COFF = 10;

		string catFile = GlassCat::CatPath + cat.name + ".AGF";
		string content = Common::ReadToString(catFile);
		istringstream ifs(content);

		cat.vecGlass.clear();
		
		while (std::getline(ifs, line))
		{
			if (InList(line, "CC")) {
				istrstream icmt(line.c_str());
				icmt >> mnemonic >> catComment;
			}
			else if (InList(line,"NM")) {
				std::string lines[6];
				lines[0] = line;

				int istart = 1; 
				std::string linex;
				std::getline(ifs, linex);
				
				if (InList(linex, "GC")) {
					istart = 1;
					istrstream icmt(linex.c_str());
					icmt >> mnemonic >> gComment;
				}
				else {
					istart = 2;
					lines[1] = linex;
				}
				
				for (int i = istart; i < 6; i++) {
					std::getline(ifs, lines[i]);
				}

				istrstream line0(lines[0].c_str());
				istrstream line1(lines[1].c_str());
				istrstream line2(lines[2].c_str());
				istrstream line3(lines[3].c_str());
				istrstream line4(lines[4].c_str());
				istrstream line5(lines[5].c_str());

				line0 >> mnemonic >> Name >> Formula >> MILID >> nd >> vd >> excsub >> status >> meltfreq;
				line1 >> mnemonic >> TCE >> TCE2 >> density >> dPgF >> IngExp;
				line2 >> mnemonic >> C[0] >> C[1] >> C[2] >> C[3] >> C[4] >> C[5] >> C[6] >> C[7] >> C[8] >> C[9];
				line3 >> mnemonic >> D0 >> D1 >> D2 >> E0 >> E1 >> Ltk >> TempUsed;
				line4 >> mnemonic >> RelCost >> CR >> FR >> SR >> AR >> PR;
				line5 >> mnemonic >> lmin >> lmax;

				Glass g(Name, cat.name);
				g.Formula = Formula;
				g.nd = nd;
				g.vd = vd;
				g.ExcSub = excsub;
				g.status = status;
				g.meltfreq = meltfreq;
				g.comment = gComment;
				g.TCE = TCE;
				g.TCE2 = TCE2;
				for (int i = 0; i < N_COFF; i++) {
					g.C[i] = C[i];
				}
				g.D0 = D0;
				g.D1 = D1;
				g.D2 = D2;
				g.E0 = E0;
				g.E1 = E1;
				g.Ltk = Ltk;
				g.Temp = TempUsed;

				g.RelCost = RelCost;
				g.CR = CR;
				g.FR = FR;
				g.SR = SR;
				g.AR = AR;
				g.PR = PR;

				g.lmin = lmin;
				g.lmax = lmax;

				cat.vecGlass.push_back(g);
			}
		}
		cat.number = cat.vecGlass.size();
	}

	void GlassCat::ExportCSV(string path) {
		string csvFile = path + name + ".csv";
		ofstream out(csvFile.c_str());
		const vector<string> params = { "name","formula","mil#","nd","vd","excsub","status","meltfreq","comments","TCE","TCE2","density","dPgF","ingExp","C0","C1","C2","C3","C4","C5","C6","C7","C8","C9","D0","D1","D2","E0","E1","Ltk","TempUsed","RelCost","CR","FR","SR","AR","PR","lambda1","lambda2" };

		string oline = "";
		for (int i = 0; i < params.size(); i++) {
			oline += params[i] + ",";
		}
		out << oline << endl;

		for (int i = 0; i < vecGlass.size(); i++) {
			Glass g = vecGlass[i];
			out << g.name << "," << g.Formula << "," << g.milnum << "," << g.nd << "," << g.vd << "," << g.ExcSub << "," << g.status << "," << g.meltfreq << ",";
			out << comment << ",";
			out << g.TCE << "," << g.TCE2 << "," << g.density << "," << g.xPgF << "," << g.IngExp << ",";
			out << g.C[0] << "," << g.C[1] << "," << g.C[2] << "," << g.C[3] << "," << g.C[4] << "," << g.C[5] << "," << g.C[6] << "," << g.C[7] << "," << g.C[8] << "," << g.C[9] << ",";
			out << g.D0 << "," << g.D1 << "," << g.D2 << "," << g.E0 << "," << g.E1 << "," << g.Ltk << "," << g.Temp << ",";
			out << g.RelCost << "," << g.CR << "," << g.FR << "," << g.SR << "," << g.AR << "," << g.PR << ",";
			out << g.lmin << "," << g.lmax << endl;
		}
	}

	void GlassCat::ExportJson(string path) {
		string jsonFile = path + name + ".json";
		ofstream out(jsonFile.c_str());
	}
}