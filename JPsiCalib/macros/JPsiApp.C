// C++ includes
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include "finalJPsiAnalysis.cc"
#include "finalJPsiAnalysisEle.cc"
#include "invMassStudy.cc"

int main(int argc, char* argv[]) {

  int applic = 1;

  char inputFileName[150];
  if ( argc < 3 ){
    std::cout << "missing argument: " << std::endl;
    std::cout << "insert A) inputFile with list of root files" << std::endl; 
    std::cout << "insert B) sample = " << std::endl;
    std::cout << "1) signal j/psi prompt"      << std::endl;
    std::cout << "2) signal j/psi non prompt"  << std::endl;
    std::cout << "3) QCD bc->e, 20-30"         << std::endl;
    std::cout << "4) QCD bc->e, 30-80"         << std::endl;
    std::cout << "5) QCD bc->e, 80-170"        << std::endl;
    std::cout << "6) QCD em. enriched, 20-30"  << std::endl;
    std::cout << "7) QCD em. enriched, 30-80"  << std::endl;
    std::cout << "8) QCD em. enriched, 80-170" << std::endl;
    return 1;
  }
  strcpy(inputFileName,argv[1]);
  int theSample = atoi(argv[2]);

  TChain *theChain = new TChain("T1");
  char Buffer[500];
  char MyRootFile[2000];
  std::cout << "input: " << inputFileName << std::endl;
  ifstream *inputFile = new ifstream(inputFileName);

  int nfiles=1;
  while( !(inputFile->eof()) ){
    inputFile->getline(Buffer,500);
    if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))) { 
      sscanf(Buffer,"%s",MyRootFile);
      theChain->Add(MyRootFile);
      std::cout << "chaining " << MyRootFile << std::endl;
      nfiles++;
    }
  }
  inputFile->close();
  delete inputFile;
  
  if (applic==0) {
    finalJPsiAnalysis jpsiStudy(theChain);
    jpsiStudy.Loop();
  }  

  if (applic==1) {
    finalJPsiAnalysisEle jpsiStudy(theChain);
    jpsiStudy.Loop(theSample);
  }  

  if (applic==2) {
    invMassStudy jpsiStudy(theChain);
    jpsiStudy.Loop();
  }  

  return 0;

}
