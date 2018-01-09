#ifndef EP2003_84_FIG1L_H
#define EP2003_84_FIG1L_H

#include <vector>
#include <iostream>

class ep2003_084_Fig1L{
 public:
 
  //using digitizer
  static const Int_t nPoints = 21;
  Double_t xVals[nPoints] = {0.014779116465863329,
			     0.024176706827309158,
			     0.03461847389558223,
			     0.04506024096385533,
			     0.054457831325301104,
			     0.06489959839357418,
			     0.07534136546184733,
			     0.08473895582329302,
			     0.09518072289156612,
			     0.10979919678714847,
			     0.12963855421686735,
			     0.1505220883534136,
			     0.17036144578313242,
			     0.1902008032128514,
			     0.22465863453815252,
			     0.2747791164658635,
			     0.3499598393574296,
			     0.45020080321285116,
			     0.5504417670682729,
			     0.6496385542168672,
			     0.749879518072289};
  
  Double_t yVals[nPoints] = {522.561874933373,
			     301.03434256927454,
			     199.03588600014234,
			     147.6121279351077,
			     106.99192136266706,
			     86.98450383896477,
			     69.11267970278135,
			     56.19038265054258,
			     45.68281878789116,
			     37.135594721899096,
			     26.90827108015079,
			     20.413523643495008,
			     14.79153981547273,
			     11.482434709596436,
			     8.127694741337956,
			     4.570354185842621,
			     2.627767333911709,
			     1.1726753380002704,
			     0.32306682590227054,
			     0.15805085512606515,
			     0.07053224901775783};
  
  Double_t bins[nPoints+1] = {0.009558232931726862,
			      0.019999999999999934,
			      0.03044176706827298, 
			      0.03983935742971878, 
			      0.05028112449799188, 
			      0.06072289156626495, 
			      0.07012048192771073, 
			      0.08056224899598383, 
			      0.0899598393574296, 
			      0.09935742971887532,
			      0.12024096385542152,
			      0.1400803212851405, 
			      0.15991967871485888,
			      0.1797590361445783, 
			      0.19959839357429718,
			      0.24971887550200791,
			      0.2998393574297188, 
			      0.4000803212851404, 
			      0.5003212851405622, 
			      0.6005622489959838, 
			      0.6997590361445782, 
			      0.7999999999999998};
  
  ep2003_084_Fig1L(){}
  ~ep2003_084_Fig1L(){}

  Int_t getNPoints(){return nPoints;}
  Bool_t checkInVal(Int_t inVal);
  Double_t getPointXVal(Int_t inVal);
  Double_t getPointYVal(Int_t inVal);
  std::vector<Double_t> getBins();
};


Bool_t ep2003_084_Fig1L::checkInVal(Int_t inVal)
{
  if(inVal < 0 || inVal >= nPoints){
    std::cout << "requested point \'" << inVal << "\' is outside accepted range 0-" << nPoints-1 << ". return -1" << std::endl;
    return false;
  }
  return true;
}


Double_t ep2003_084_Fig1L::getPointXVal(Int_t inVal)
{
  if(!checkInVal(inVal)) return -1;
  else return xVals[inVal];
}


Double_t ep2003_084_Fig1L::getPointYVal(Int_t inVal)
{
  if(!checkInVal(inVal)) return -1;
  else return yVals[inVal];
}


std::vector<Double_t> ep2003_084_Fig1L::getBins()
{
  std::vector<Double_t> binsV;
  for(Int_t i = 0; i < nPoints+1; ++i){binsV.push_back(bins[i]);}
  return binsV;
}

#endif
