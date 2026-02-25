void CreateTree(TString particle = "B")
{
	TChain *t = new TChain("Jets/DecayTree");
	if (particle == "B")
	{
	  t->Add("Bjet_b_2016_MD_09082022_url333_0_800.root");
	  t->Add("Bjet_b_2016_MD_09082022_url333_800_1429.root");
	  t->Add("Bjet_b_2016_MU_09082022_url334_0_800.root");	  
	  t->Add("Bjet_b_2016_MU_09082022_url334_800_1299.root");
	  t->Add("Bjet_b_2017_MD_09082022_url335_0_800.root");
	  t->Add("Bjet_b_2017_MD_09082022_url335_800_882.root");
	  t->Add("Bjet_b_2017_MU_09082022_url336_0_800.root");
	  t->Add("Bjet_b_2017_MU_09082022_url336_800_846.root");
	  t->Add("Bjet_b_2018_MD_09082022_url337_0_800.root");
	  t->Add("Bjet_b_2018_MD_09082022_url337_800_988.root");
	  t->Add("Bjet_b_2018_MU_09082022_url338_0_800.root");
	  t->Add("Bjet_b_2018_MU_09082022_url338_800_1042.root");				
	}
	
	if (particle == "D")
	{
	  t->Add("Djet_2016_MD_09082022_url327_0_535.root");
	  t->Add("Djet_2016_MU_09082022_url328_0_452.root");
	  t->Add("Djet_2017_MD_09082022_url329_0_345.root");	  
	  t->Add("Djet_2017_MU_09082022_url330_0_327.root");
	  t->Add("Djet_2018_MD_09082022_url331_0_372.root");
	  t->Add("Djet_2018_MU_09082022_url332_0_392.root");			
	}	
	
	TString className = particle;
	className += "JetTree";
	t->MakeClass(className);

}
