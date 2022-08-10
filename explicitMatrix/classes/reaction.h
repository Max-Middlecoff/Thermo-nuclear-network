
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>

using namespace std;
using std::string;

/* Class Reaction to describe all reactions in the network.  Create new instance
 * of this class for each reaction in the network. Inherits from Utilities */

#define THIRD 0.333333333333333

class Reaction: public Utilities {
    
    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
    private:
        
        int reacIndex;               // Index of reaction in current run for each (Z,N)
        int reacClass;               // Reaction class for reaclib (1-8)
        int reacGroupClass;          // Reaction group class (1-5 surrogate for A-E)
        int rgindex;                 // Index of RG containing reaction (0, 1, ... #RG)
        int RGmemberIndex;           // Index of reaction within its reaction group
        
        char reacGroupClassLett;     // Letter equivalent (A-E) for reacGroupClass
        char* reacGroupSymbol;       // Schematic equil reaction (e.g. a+b<->c)
        int numberReactants;         // Number species on the left side of reaction
        int numberProducts;          // Number species on the right side of reaction
        int numberIsotopes;          // numberReactants + numberProducts in reaction
        char* reacString;            // String describing reaction
        char resonanceType;          // Whether resonant (r) or non-resonant (nr)
        int isEC;                    // Whether electron capture reaction (1) or not (0)
        int isReverse;               // Whether reverse reaction (1) or not (0)
        int ispeforward;             // Whether reactions is labeled "forward" in PE scheme
        bool isEquil;                // Whether in a RG satisfying PE conditions
        
        double Q;                    // Q-value for reaction
        double prefac;               // The eta prefac for rates
        double p[7];                 // ReacLib parameters
        int reactantZ[3];            // Array holding Z of reactants
        int reactantN[3];            // Array holding N of reactants
        int reactantA[3];            // A = Z+N
        int productZ[4];             // Array holding Z of products
        int productN[4];             // Array holding N of products
        int productA[4];             // A = Z+N
        int reactantIndex[3];        // Index of species isotope vector for each reactant isotope
        int productIndex[4];         // Index of species isotope vector for each product isotope
        int isoIndex[7];             // Index of species isotope vector for all isotopes in reaction
        
        // Precomputed temperature factors for ReacLib rates.  Computed in computeTfacs(T9), where
        // T9 is the temperature in units of 10^9 K.
        
        double T93;                  // T9^{1/3}
        double t1;                   // 1/T9
        double t2;                   // 1/T9^{1/3}
        double t3;                   // T93
        double t4;                   // T9
        double t5;                   // T9^{5/3} = (T93)^5
        double t6;                   // ln (T9)
        
        double Dens[3];              // Powers of density
        double densfac;              // prefac x powers of density
        double rate;                 // Current T-dependent part of rate for reaction
        double Rrate;                // Current full rate (s^-1) for reaction
        double flux;                 // Current effective net flux of reaction if PE
        double flux_true;            // Current true net flux of reaction
        double dErate;               // Current rate of energy release
        char cs[20];                 // Utility character array
        char ccs[20];                // Utility character array
        char* ss;                    // Utility string
  
  
    public:
        
        // Constructor executed when objects are instantiated
        
        Reaction(bool reacIsActive[]){
            
            // Set all reaction objects to not-equilibrated initially.
            
            isEquil = false;
            reacIsActive[reacIndex] =  true;

        }
        
        // Public Reaction setter functions to set values of private class fields
        
        void setreacIndex(int index){reacIndex = index; }
        
        void setreacClass(int index){reacClass = index; }
        
        void setreacGroupClass(int index, int RGclass[], string RGstring[]){
            
            reacGroupClass = index;         // Field in this class
            RGclass[reacIndex] =  index;    // In main
            
            // Use this index to set letters in reacGroupClassLett
            // and symbols in reacGroupSymbol
            
            switch(index){
                case 1:
                    reacGroupClassLett = 'A';
                    reacGroupSymbol = Utilities::stringToChar("a<->b");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 2:
                    reacGroupClassLett = 'B';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 3:
                    reacGroupClassLett = 'C';
                    reacGroupSymbol = Utilities::stringToChar("a+b+c<->d");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 4:
                    reacGroupClassLett = 'D';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c+d");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
                case 5:
                    reacGroupClassLett = 'E';
                    reacGroupSymbol = Utilities::stringToChar("a+b<->c+d+e");
                    RGstring[reacIndex] = reacGroupSymbol;
                    break;
            }
        }
        
        void setreacGroupSymbol(char* s){reacGroupSymbol = s; }
        
        void setRGmemberIndex(int i, int RGMemberIndex[]){
            
            RGmemberIndex = i;                   // Field in this class
            RGMemberIndex[reacIndex] = i;        // Write to corresponding array 
            
        }
        
        void setrgindex(int i){rgindex = i; }
        
        void setnumberReactants(int i, int NumReactingSpecies[]){
            
            numberReactants = i;                 // Field in this class
            NumReactingSpecies[reacIndex] = i;   // Array
            
        }
        
        void setnumberProducts(int i, int NumProducts[]){
            
            numberProducts = i;                  // Field in this class
            NumProducts[reacIndex] = i;          // Array
            
        }
        
        void setreacString(
            char* s,
            char reacLabel[][35] // @todo: SHOULD NOT DO THIS!!!!!
        ){ 
            
            // Set field of this class
            
            reacString = s; 
            
            // Set corresponding character array reacLabel 
            
            char p[20];  
            for (int i = 0; i < sizeof(p); i++) { 
                p[i] = s[i]; 
                reacLabel[reacIndex][i] = p[i];
            }
        }
        
        void setisEC(int i){ isEC = i; }
        
        void setisReverse(int i){ isReverse = i; }
        
        void setispeforward(int i, int isPEforward[]){
            ispeforward = i;
            isPEforward[reacIndex] = i;
        }
        
        void setisEquil(bool e){isEquil = e;}
        
        void setQ(double q){ Q = q; }
        
        void setprefac(double p){ prefac = p; }
        
        void setp(double P[]){
            for(int i=0; i<7; i++){
                p[i] = P[i];
            }
        }
        
        void setreactantZ(int z[], int reacZ[][4] ){
            for (int i=0; i<numberReactants; i++){
                
                // Change field in this class (Reaction)
                
                reactantZ[i] = z[i];
                
                // Change array 
                
                reacZ[getreacIndex()][i] = z[i];
            }

        }
        
        void setreactantN(int n[], int reacN[][4] ){
            for (int i=0; i<numberReactants; i++){
                
                // Change field in this class (Reaction)
                
                reactantN[i] = n[i];
                
                // Change array
                
                reacN[getreacIndex()][i] = n[i];
            }
        }
        
        void setreactantA(void){
            for (int i=0; i<numberReactants; i++){
                reactantA[i] = reactantZ[i] + reactantN[i];
            }
        }
        
        void setproductZ(int z[], int prodZ[][4] ){
            for (int i=0; i<numberProducts; i++){
                
                // Change field in this class (Reaction)
                
                productZ[i] = z[i];
                
                // Change array
                
                prodZ[getreacIndex()][i] = z[i];
            }
        }
        
        void setproductN(int n[], int prodN[][4] ){
            for (int i=0; i<numberProducts; i++){
                
                // Change field in this class (Reaction)
                
                productN[i] = n[i];
                
                // Change array
                
                prodN[getreacIndex()][i] = n[i];
            }
        }
        
        void setproductA(void){
            for (int i=0; i<numberProducts; i++){
                productA[i] = productZ[i] + productN[i];
            }
        }
        
        void setreactantIndex(int n[], int ReactantIndex[][4] ){
            for (int i=0; i<numberReactants; i++){
                reactantIndex[i] = n[i];              // Class field
                ReactantIndex[reacIndex][i] = n[i];   // Array
            }
        }
        
        void setproductIndex(int n[], int ProductIndex[][4] ){
            for (int i=0; i<numberProducts; i++){
                productIndex[i] = n[i];               // Class field
                ProductIndex[reacIndex][i] = n[i];    // Array
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes no arguments
        // and constructs isoIndex[] as the concatenation of reactantIndex[]
        // and productIndex[], assuming that those fields have been populated.
        
        void setisoIndex(void){
            for(int i=0; i<numberReactants; i++){
                isoIndex[i] = reactantIndex[i];
            }
            for(int i=0; i<numberProducts; i++){
                isoIndex[i+numberReactants] = productIndex[i];
            }
        }
        
        // Overloaded versions of setisoIndex.  This version takes two 
        // arguments and sets a value of a particular isoIndex[].
        
        void setisoIndex(int i, int j){isoIndex[i] =  j;}
        
        void setdensfac(double d){ densfac = d;}
        
        void setrate(double r){ rate = r; }
        
        void setRrate(double r){ Rrate = r; }
        
        void setflux(double f){ flux = f; }
        
        void setflux_true(double f){ flux_true = f;}
        
        // Overloaded dErate setter
        
        void setdErate(double r){ dErate = r; }
        void setdErate(void){ dErate = flux_true*Q; }
        
        
        // Public Reaction getter functions to get values in private fields
        
        int getreacIndex(){ return reacIndex; }

        int getreacClass(){ return reacClass; }
        
        int getreacGroupClass(){ return reacGroupClass; }
        
        // Return reacGroupSymbol as string
        
        char* getreacGroupSymbol(){ return reacGroupSymbol; }
        
        int getRGmemberIndex(){return RGmemberIndex;}
        
        int getrgindex(){return rgindex;}
        
        char getreacGroupClassLett(){ return reacGroupClassLett; }
        
        int getnumberReactants(){ return numberReactants; }
        
        int getnumberProducts(){ return numberProducts; }
        
        // return reacString as string
        
        char* getreacString(){ return reacString; }
        
        // The function getreacChar() returns the string reacString as a
        // pointer to a character array that will work in printf. Alternatively,
        // Utilities::stringToChar() will do same thing.
        
        char* getreacChar(){
            return ccs;
        }
        
        int getisEC(){ return isEC; }
        
        int getisReverse(){ return isReverse; }
        
        int getispeforward(){return ispeforward;}
        
        bool getisEquil(){return isEquil;}
        
        double getQ(){ return Q; }
        
        double getprefac(){ return prefac; }
        
        double getp(int k){ return p[k]; }
        
        int getreactantZ(int k){
            
            if(k > numberReactants-1){
                printf("\n\nERROR: k=%d larger than numberReactants-1 %d; Stopping\n", 
                    k, numberReactants);
                exit(-1);                   // Exit program since something is corrupt
            } else {
                return reactantZ[k];
            }
        }
        
        int getreactantN(int k){
            if(k > numberReactants-1){
                printf("\n\nERROR: k-1=%d larger than number reactants %d", 
                    k, numberReactants);
                return -1;
            } else {
                return reactantN[k];
            }
        }
        
        int getreactantA(int i){
           return (reactantZ[i] + reactantN[i]);
        }
        
        int getproductZ(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k=%d larger than number products %d", 
                    k, numberProducts);
                return -1;
            } else {
                return productZ[k];
            }
        }
        
        int getproductN(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d", 
                    k, numberProducts);
                return -1;
            } else {
                return productN[k];
            }
        }
        
        int getproductA(int i){
            return (productZ[i] + productN[i]);
        }
        
        int getreactantIndex(
                                int k,
                                char reacLabel[][35] // @todo: SHOULD NOT DO THIS!!!!!
                            ){
            if(k > numberReactants-1){
                printf("\nReac=%d %s", reacIndex, reacLabel[reacIndex]);
                ss = Utilities::stringToChar(
                    "\nERROR Reaction::getreactantIndex(k): k = %d larger than #reactants-1 = %d");
                printf(stringToChar(ss), k, numberReactants-1);
                return -1;
            } else {
                return reactantIndex[k];
            }
        }
        
        int getproductIndex(int k){
            if(k > numberProducts-1){
                printf("\n\nERROR: k-1=%d larger than number products %d", 
                    k, numberProducts);
                return -1;
            } else {
                return productIndex[k];
            }
        }
        
        double getdensfac(){ return densfac; }
        
        int getisoIndex(int i){return isoIndex[i];}
        
        double getrate(){ return rate; }
        
        double getRrate(){ return Rrate; }
        
        double getflux(){ return flux; }
        
        double getflux_true(){ return flux_true;}
        
        double getdErate(){ return dErate; }
        
        
        // Function Reaction::setupFplusFminus() to set up F+ and F- index for each
        // isotope and to find non-vanishing F+ and F- source terms in network.
        // Function declared static so it can be called as Reaction::setupFplusFminus() 
        // without having to instantiate.
        
        static void setupFplusFminus(
                                        int FplusIsotopeCut[], 
                                        int numFluxPlus[], 
                                        int FminusIsotopeCut[], 
                                        int numFluxMinus[],
                                        int numberSpecies,
                                        int totalFplus,
                                        int FplusIsotopeIndex[],
                                        int FminusIsotopeIndex[],
                                        int totalFminus,
                                        FILE * pFileD,
                                        int MapFplus[],
                                        int MapFminus[],
                                        int tempInt1[],
                                        int tempInt2[],
                                        char reacLabel[][35], // @todo: SHOULD NOT DO THIS!!!!!
                                        int FplusMin[],
                                        int FplusMax[],
                                        int FminusMin[],
                                        int FminusMax[],
                                        int ISOTOPES,
                                        int SIZE,
                                        int reacMask[][48], // @todo: SHOULD NOT DO THIS!!!!!
                                        double FplusFac[],
                                        double FminusFac[]
                                    ){
           
            FplusIsotopeCut[0] = numFluxPlus[0];
            FminusIsotopeCut[0] = numFluxMinus[0];
            for(int i=1; i<numberSpecies; i++){
                FplusIsotopeCut[i] = numFluxPlus[i] + FplusIsotopeCut[i-1];
                FminusIsotopeCut[i] = numFluxMinus[i] + FminusIsotopeCut[i-1];
            }
            
            int currentIso = 0;
            for(int i=0; i<totalFplus; i++){
                FplusIsotopeIndex[i] = currentIso;
                if(i == (FplusIsotopeCut[currentIso]-1)) currentIso ++;
            }
            
            currentIso = 0;
            for(int i=0; i<totalFminus; i++){
                FminusIsotopeIndex[i] = currentIso;
                if(i == (FminusIsotopeCut[currentIso]-1)) currentIso ++;
            }
            
            fprintf(pFileD,"\n\nMapFplus:\n");
            
            for(int i=0; i<totalFplus; i++){
                MapFplus[i] = tempInt1[i];
                fprintf(pFileD, "\ni=%d MapFplus[i]=%d  %s",
                    i, MapFplus[i], reacLabel[MapFplus[i]]);
            }
            
            fprintf(pFileD,"\n\nMapFminus:\n");
            
            for(int i=0; i<totalFminus; i++){
                MapFminus[i] = tempInt2[i];
                fprintf(pFileD, "\ni=%d MapFminus[i]=%d  %s",
                        i, MapFminus[i], reacLabel[MapFminus[i]]);
            }
            
            // Populate the FplusMin and FplusMax arrays
            
            FplusMin[0] = 0;
            FplusMax[0] = numFluxPlus[0]-1;
            
            for(int i=1; i<numberSpecies; i++){
                FplusMin[i] = FplusMax[i-1] + 1;
                FplusMax[i] = FplusMin[i] + numFluxPlus[i] -1 ;	
            }
            
            // Populate the FminusMin and FminusMax arrays
            
            FminusMin[0] = 0;
            FminusMax[0] = numFluxMinus[0]-1;
            for(int i=1; i<numberSpecies; i++){
                FminusMin[i] = FminusMax[i-1] + 1;
                FminusMax[i] = FminusMin[i] + numFluxMinus[i] -1 ;	
            }
            
            // Populate the FplusFac and FminusFac arrays that hold the factors counting the
            // number of occurrences of the species in the reaction.  Note that this can only
            // be done after ReactionVector::parseF() has been run to give reacMask[i][j].
            
            int tempCountPlus = 0;
            int tempCountMinus = 0;
            for(int i=0; i<ISOTOPES; i++){
                for(int j=0; j<SIZE; j++)
                {
                    if(reacMask[i][j] > 0)
                    {
                        FplusFac[tempCountPlus] = (double)reacMask[i][j];
                        tempCountPlus ++;
                    }
                    else if(reacMask[i][j] < 0)
                    {
                        FminusFac[tempCountMinus] = -(double) reacMask[i][j];
                        tempCountMinus ++;
                    }	
                }
            }
        }
        
        
        // Static function Reaction::populateFplusFminus() to populate F+ and F- for each
        // isotope set up in setupFplusFminus() from master flux array. Function declared
        // static so it can be called as Reaction::populateFplusFminus() without having
        // to instantiate.
        
        static void populateFplusFminus(
                                         int totalFplus,
                                         int MapFplus[],
                                         double Fplus[],
                                         double FplusFac[],
                                         int totalFminus,
                                         int MapFminus[],
                                         double Flux[],
                                         double Fminus[],
                                         double FminusFac[]
            ){
            
            // Populate the F+ and F- arrays from the master Flux array
            
            for(int i=0; i<totalFplus; i++){
                int indy = MapFplus[i];
                Fplus[i] = FplusFac[i]*Flux[indy];
            }
            
            for(int i=0; i<totalFminus; i++){
                int indy = MapFminus[i];
                Fminus[i] = FminusFac[i]*Flux[indy];
            }
        }
        
        
        // Reaction::computeConstantFacs() to compute the reaction factors 
        // that are constant for a given temperature and density.  Reset each 
        // time T9 or rho changes but stays constant as long as they don't change.  
        // This is required at the beginning of each network integration of the 
        // hydro timestep, since the density and temperature (in addition to
        // abundances because of advection) will generally change over a hydro timestep 
        // in each zone.
        
        void computeConstantFacs(double T9, double rho){

            // Temperature factors in ReacLib rate formula.

            T93 = powf(T9, THIRD); 
            t1 = 1/T9;
            t2 = 1/T93;
            t3 = T93;
            t4 = T9;
            t5 = T93*T93*T93*T93*T93;
            t6 = logf(T9);
            
            // Multiply the statistical prefactor by the appropriate 
            // density factors (1 for 1-body, rho for 2-body, and rho^2 
            // for 3-body reactions).
            
            Dens[0] = 1.0f;
            Dens[1] = rho;
            Dens[2] = rho*rho;
            densfac = prefac * Dens[numberReactants - 1];
            setdensfac(densfac); 
            
        }
        
        // Function Reaction::computeTfacs(double) to set temperature factors in ReacLib
        // rate formula.  Reset each time T9 changes but
        // stays constant as long as T9 doesn't change.
        
        void computeTfacs(double T9){
            T93 = powf(T9, THIRD); 
            t1 = 1/T9;
            t2 = 1/T93;
            t3 = T93;
            t4 = T9;
            t5 = T93*T93*T93*T93*T93;
            t6 = logf(T9);
        }
        
        // Reaction::computeDensityFactors(double) to multiply the statistical prefactor by 
        // the appropriate density factors (1 for 1-body, rho for 2-body, and rho^2 for 3-body. 
        // This is required at the beginning of each network integration of the hydro timestep, 
        // since the density will generally change over a hydro timestep in each zone.
        
        void computeDensityFactors(double rho){
            
            Dens[0] = 1.0f;
            Dens[1] = rho;
            Dens[2] = rho*rho;
            densfac = prefac * Dens[numberReactants - 1];
            setdensfac(densfac);
        }
        
        // Reaction::computeRate(double, double) to compute rates at T and rho. 
        // The quantity rate is the temperature-dependent part, including a possible 
        // partition-function correction.  The quantity Rrrate is rate multiplied by
        // appropriate density and statistical factors, which give units of s^-1.  The 
        // flux follows from multiplying Rrate by appropriate abundances Y in computeFlux().
        
        void computeRate(double T9, double rho, double Rate[], bool dopf, double pfCut9){
            
            // Temperature-dependent rate from ReacLib library
            
            rate = expf( p[0] + t1*p[1] + t2*p[2] + t3*p[3] + t4*p[4] + t5*p[5] + t6*p[6] );
            
            // Correct rate by multiplying by partition functions and store
            
            pfUpdate(dopf, T9, pfCut9);
            setrate(rate);

            // Full rate factor in s^-1 (rate above multiplied by density factors)
            
            Rrate = getdensfac() * rate;
            setRrate(Rrate);
            
            // Write to rate array
            
            Rate[getreacIndex()] = Rrate;
            
        }
        
        
        // Function Reaction::pfUpdate() to correct the rates using partition 
        // function factors if appropriate.
        
        void pfUpdate(bool dopf, double T9, double pfCut9){
            
            double pfnum;
            double pfden;
            double pfFactor;
            
            // Make a partition function correction if this is reverse reaction in
            // sense defined in ReacLib (defined by field Reaction::isReverse=true). 
            // Realistic calculations at higher temperatures should use
            // the partition functions so generally set dopf=true unless testing.
            // Partition functions are very near 1.000 if T9 < 1, so we will typically
            // only implement partition function correction if T9 > pfCut9 = 1.0, but
            // the table of partition functions allows pfCut9 as small as 0.1.
            // Interpolation is in the log10 of the temperature, so pass log10(T9)
            // rather than T9 to pfInterpolator (index, logt9). Because for the 
            // temperatures of interest the partition functions for all light ions
            // (protons, neutrons, alphas, tritons) are equal to 1.0, the structure
            // of the 8 Reaclib reaction classes specified by Reaction::reacClass
            // means that this correction is only required for reverse reactions
            // in Reaclib classes reacClass = 2, 5.
            
            if(dopf && T9 > pfCut9 && isReverse){
                
                if(reacClass == 2){
                    pfden = pfInterpolator (reactantIndex[0], log10(T9));
                    pfnum = pfInterpolator (productIndex[1], log10(T9));
                } else if(reacClass == 5){
                    pfden = pfInterpolator (reactantIndex[1], log10(T9));
                    pfnum = pfInterpolator (productIndex[1], log10(T9));
                } else {
                    pfden = 1.0;
                    pfnum = 1.0;
                }
                pfFactor = pfnum/pfden;
                rate *= pfFactor;
                
            }
        }
        
        
        // ------------------------------------------------------------------------
        // Reaction::pfInterpolator(int, double) to return
        // partition function of isotope labeled by isoIndex at log_10 of
        // temperature T9. Note that the 2nd argument is log10(T9), not T9,
        // because the interpolation in the partition function table is in the 
        // log10 of the temperature.  The following commented-out code assumes
        // that the object interpolatepf of the SplineInterpolator class has
        // first invoked the interpolatepf.bisection method to use bisection 
        // to find the interval containing root and store the lower index of
        // that interval in lowPFindex. Then SplineInterpolator interpolates
        // the root restricted to that interval.  This guards against the
        // spline interpolator finding the wrong root if there are multiple
        // roots (as could be true in the general case, though probably not here
        // since the function is typically monotonic).
        // ------------------------------------------------------------------------
        
        double pfInterpolator(int index, double logt9) {
            
            // Following commented out for testing purposes until spline interpolator
            // is implemented
            
//             double rdt;
//             double term1;
//             double term2;
//             double sumterms;
//             double bob;
//             rdt = (logt9 - Tpf[lowPFindex]) / (Tpf[lowPFindex + 1] - Tpf[lowPFindex]);
//             term1 = rdt * Math.log(pf[Z][N][lowPFindex + 1]);
//             term2 = (1.0 - rdt) * Math.log(pf[Z][N][lowPFindex]);
//             sumterms = term1 + term2;
//             bob = Math.exp(sumterms);
//             // System.out.println("PF stuff: "+t9+" "+Z+" "+N+" "+rdt+" "+sumterms+" "+bob);
//             return bob;
            
            return 1.0;  // Temporary
            
        }
        
        // Function Reaction::showRates() to display computed rates for this
        // Reaction object.
        
        void showRates( FILE * pFileD ){
            fprintf(pFileD, "\n%d %19s RG=%d densfac=%6.3e rate= %8.5e Rrate=%8.5e", 
                getreacIndex(), getreacChar(), getreacGroupClass(), getdensfac(), 
                getrate(), getRrate()
            );
        }
        
        
        // Function Reaction::computeFlux() to compute fluxes for reactions.  This is
        // where net flux for a reaction group is set to zero if reaction group is 
        // in equilibrium. The computed flux will be set in the flux field of the
        // Reaction object, and also in the array Flux[].
        
        void computeFlux( bool reacIsActive[], double Flux[], double Y[], 
                            double fastestCurrentRate, 
                            int fastestCurrentRateIndex, 
                            double slowestCurrentRate, 
                            int slowestCurrentRateIndex,
                            double fastestOverallRate, 
                            int fastestOverallRateIndex, 
                            double slowestOverallRate,
                            int slowestOverallRateIndex,
                            double timeMaxRate,
                            double t   
                        ){
            
            // If this reaction is in a RG that is not in equilibrium, we
            // need to compute its flux.  The formula for the flux will depend on
            // whether the reaction is 1-body, 2-body, or 3-body, and is selected by
            // the switch statement in fluxChooser(numberReactants).
            
            isEquil = !reacIsActive[reacIndex];  // Set isEquil in this Reaction object
            
            if( reacIsActive[reacIndex] ) {
                
                fluxChooser(numberReactants, Y, Flux, 
                            fastestCurrentRate, 
                            fastestCurrentRateIndex, 
                            slowestCurrentRate, 
                            slowestCurrentRateIndex,
                            fastestOverallRate, 
                            fastestOverallRateIndex, 
                            slowestOverallRate,
                            slowestOverallRateIndex,
                            timeMaxRate,
                            t );
            
            // Otherwise the reaction is in an equilibrated reaction group, so set its
            // flux identically to zero.
            
            } else {
                
                fluxChooser(numberReactants, Y, Flux, 
                            fastestCurrentRate, 
                            fastestCurrentRateIndex, 
                            slowestCurrentRate, 
                            slowestCurrentRateIndex,
                            fastestOverallRate, 
                            fastestOverallRateIndex, 
                            slowestOverallRate,
                            slowestOverallRateIndex,
                            timeMaxRate,
                            t );
                flux = 0.0;                 // Put in present (Reaction) object flux field
                Flux[reacIndex] = 0.0;      // Put in main flux array

            }
            
        }   // End of function computeFlux()

  
  
    // Function Reaction::fluxChooser(int) to switch among flux formulas for
    //  one-body, two-body, or three-body reactions.
      
    void fluxChooser(
                            int bodies, 
                            double Y[], 
                            double Flux[], 
                            double fastestCurrentRate, 
                            int fastestCurrentRateIndex, 
                            double slowestCurrentRate, 
                            int slowestCurrentRateIndex,
                            double fastestOverallRate, 
                            int fastestOverallRateIndex, 
                            double slowestOverallRate,
                            int slowestOverallRateIndex,
                            double timeMaxRate,
                            double t
                    ){

        double kfac;
            
            switch(bodies){
                
                case 1:    // 1-body reactions
                    
                    kfac = Rrate;
                    flux = kfac*Y[ reactantIndex[0] ];	 // In Reaction object flux field
                    Flux[reacIndex] = flux;              // In main flux array
                    fastSlowRates( 
                                    kfac, 
                                    fastestCurrentRate, 
                                    fastestCurrentRateIndex,  
                                    slowestCurrentRate, 
                                    slowestCurrentRateIndex, 
                                    fastestOverallRate, 
                                    fastestOverallRateIndex,  
                                    slowestOverallRate, 
                                    slowestOverallRateIndex,
                                    t,
                                    timeMaxRate
                                );
                    break;
                    
                case 2:	   // 2-body reactions
                    
                    kfac = Rrate * Y[ reactantIndex[0] ];
                    flux = kfac * Y[ reactantIndex[1] ];  // Put in Reaction object flux field
                    Flux[reacIndex] = flux;               // Put in main flux array
                    fastSlowRates(kfac, 
                                    fastestCurrentRate, 
                                    fastestCurrentRateIndex,  
                                    slowestCurrentRate, 
                                    slowestCurrentRateIndex, 
                                    fastestOverallRate, 
                                    fastestOverallRateIndex,  
                                    slowestOverallRate, 
                                    slowestOverallRateIndex,
                                    timeMaxRate,
                                    t
                                );
                    
                    break;
                    
                case 3:	   // 3-body reactions
                    
                    kfac = Rrate * Y[ reactantIndex[0] ] * Y[ reactantIndex[1] ];
                    flux = kfac * Y[ reactantIndex[2] ];  // Put in Reaction object flux field
                    Flux[reacIndex] = flux;               // Put in main flux array
                    fastSlowRates(kfac, 
                                    fastestCurrentRate, 
                                    fastestCurrentRateIndex,  
                                    slowestCurrentRate, 
                                    slowestCurrentRateIndex, 
                                    fastestOverallRate, 
                                    fastestOverallRateIndex,  
                                    slowestOverallRate, 
                                    slowestOverallRateIndex,
                                    timeMaxRate,
                                    t
                                );
                    
                    break;
                    
            }  // End switch 

    }       // end function Reaction::fluxChooser(int)
      
      
        // Reaction::fastSlowRates(double) to store fastest and slowest rates 
        // in this timestep, and overall in the calculation. These rates are the
        // rates kfac computed in computeFlux().
        
        void fastSlowRates(
                            double testRate, 
                            double fastestCurrentRate, 
                            int fastestCurrentRateIndex, 
                            double slowestCurrentRate, 
                            int slowestCurrentRateIndex,
                            double fastestOverallRate, 
                            int fastestOverallRateIndex, 
                            double slowestOverallRate,
                            int slowestOverallRateIndex,
                            double timeMaxRate,
                            double t
                          ){
            
            if (testRate > fastestCurrentRate) {
                fastestCurrentRate = testRate;
                fastestCurrentRateIndex = getreacIndex();
            }
            
            if (testRate < slowestCurrentRate && testRate > 0.0) {
                slowestCurrentRate = testRate;
                slowestCurrentRateIndex = getreacIndex();
            }
            
            if (testRate > fastestOverallRate) {
                fastestOverallRate = testRate;
                fastestOverallRateIndex = getreacIndex();
                timeMaxRate = t;
            }
            
            if (testRate < slowestOverallRate){
                slowestOverallRate = testRate;
                slowestOverallRateIndex = getreacIndex();
            }
            
        }   // End of function Reaction::fastSlowRates()
                
};  // End class Reaction

