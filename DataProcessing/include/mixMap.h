#ifndef MIXMAP
#define MIXMAP

class MixMap{

  public:
    MixMap(unsigned char maxMult = 200);
    ~MixMap();
    
    void addElement(unsigned char mult, int indx);
    int getNextElement(unsigned char mult, int indx);
    void setCurrentElement(unsigned char mult, int indx);

  private:
    void buildRanges();
    void resetCurrEvt(unsigned char mult);
    inline unsigned char multCheck(unsigned char mult);

    bool isSorted;
    unsigned char maxM;
    std::vector< std::pair< unsigned char, int >> multMap;
    std::vector< int > nEvtsPerMult;
    std::vector< std::pair< int, int > > multRange;
    std::vector< int > multCurrEvt;

};

//for multiplicities above maxMult, bin everything together so it's easy to match them
inline unsigned char MixMap::multCheck( unsigned char mult){ return (mult > maxM-1) ? maxM-1 : mult;}

void MixMap::setCurrentElement(unsigned char mult, int indx){
  mult = multCheck(mult);
 if(isSorted == false){
    std::cout << "Sorting mixing map..." << std::endl;
    std::sort(multMap.begin(), multMap.end());
    buildRanges();
    std::cout << "Done sorting map..." << std::endl;
    isSorted = true;
  }
  if(nEvtsPerMult.at(mult)<2){
    //std::cout << "Warning: less than 2 events in sample passing matching criterion.  Cannot mix!  Returning the events own index: " << indx << std::endl;
    return;
  }

  std::pair< unsigned char, int> tempPair = std::pair< unsigned char, int >(mult, indx);
  std::vector< std::pair< unsigned char, int> >::iterator start = multMap.begin()+multRange.at(mult).first;
  std::vector< std::pair< unsigned char, int> >::iterator end = multMap.begin()+multRange.at(mult).second;
  auto it = std::find(start, end, tempPair);

  if(it==end) return;

  int index = it-multMap.begin();
  multCurrEvt.at(mult) = index;

}

void MixMap::buildRanges(){
  int boundary1 = 0;
  int boundary2 = nEvtsPerMult.at(0);
  std::vector< std::pair<  int, int >  > tempVec;
  std::vector< int > tempMultCurrEvt; 

  for(int i = 0; i<maxM; i++){
    tempMultCurrEvt.push_back(boundary1);
    tempVec.push_back(std::pair< int, int>(boundary1, boundary2));
    boundary1 = boundary2;
    if(i<maxM-1) boundary2 += nEvtsPerMult.at(i+1);
  }
 
  multCurrEvt.swap(tempMultCurrEvt);
  multRange.swap(tempVec);
}

void MixMap::resetCurrEvt(unsigned char mult){
  multCurrEvt.at(mult) = multRange.at(mult).first;
}

int MixMap::getNextElement(unsigned char mult, int indx){
 mult = multCheck(mult); 
 if(isSorted == false){
    std::cout << "Sorting mixing map..." << std::endl;
    std::sort(multMap.begin(), multMap.end());
    buildRanges();
    std::cout << "Done sorting map..." << std::endl;
    isSorted = true;
  }
  if(nEvtsPerMult.at(mult)<2){
    //std::cout << "Warning: less than 2 events in sample passing matching criterion.  Cannot mix!  Returning the events own index: " << indx << std::endl;
    return indx;
  }

  int mixIndx = multMap.at(multCurrEvt.at(mult)).second;
  multCurrEvt.at(mult)++;
  if(multCurrEvt.at(mult) >= multRange.at(mult).second) resetCurrEvt(mult);
  
  if(mixIndx != indx) return mixIndx;
  else{
    mixIndx = multMap.at(multCurrEvt.at(mult)).second;
    multCurrEvt.at(mult)++;
    if(multCurrEvt.at(mult) >= multRange.at(mult).second) resetCurrEvt(mult);
    return mixIndx;
  } 
}

void MixMap::addElement(unsigned char mult, int indx){
  isSorted = false;
  mult = multCheck(mult);
  multMap.push_back(std::pair< unsigned char, int>(mult, indx));
  nEvtsPerMult.at(mult)++;
}

MixMap::MixMap(unsigned char maxMult){
  isSorted = false;
  maxM = maxMult+1;
  nEvtsPerMult = std::vector< int >(maxM, 0);
}

MixMap::~MixMap(){}


#endif
