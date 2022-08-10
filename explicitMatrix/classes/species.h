
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <ctime>
#include <array>

using namespace std;
using std::string;

/* Class Species to describe all species in the network.  Create new instance
 * for each isotope in the network. Inherits from class Utilities. */

template<std::size_t PF>
class Species: public Utilities {
    public:
        static const std::size_t size = PF;
        
        // Static method to to renorm all X to sum of one on read-in if
        // invoked after readin of network. Now we start with sumX=1 to machine precision, 
        // even if the conversion of the Ys read in to Xs does not give sumX = 1 to machine
        // precision.
        
        // Public setter functions to set values of private Species class fields
        
        void setisoIndex(int iso){isoindex = iso;}
        
        void setZ(int z, int Z[]){
            ZZ = z;             // Class field
            Z[isoindex] = z;    // Array
        }
        
        void setN(int n, int N[]){
            NN = n;             // Class field
            N[isoindex] = n;    // Array
        }
        
        void setA(int a, int AA[]){
            A = a;              // Class field
            AA[isoindex] = a;   // Array
        }
        
        void setY(double y, double Y[], double X[]){
            YY = y;             // Class field
            XX = (double)A*YY;  // Corresponding X
            Y[isoindex] = y;    // Array
            X[isoindex] = XX;   // Array
        }
        
        void setY0(double y){ YY0 = y; }
        
        void setX(double x, double X[]){
            XX = x;             // Class field
            YY = XX/(double)A;  // Corresponding Y
            X[isoindex] = x;    // Array
        }
        
        void setM(double m, double massExcess[]){
            MassExcess = m;             // Class field
            massExcess[isoindex] = m;   // Array
        }
        
        void setisoLabel(
                            char *lpt,
                            char isoLabel[][5]  // @todo: SHOULD NOT DO THIS!!!!!
                        ){
            
            // Fill character array.  *lpt points to initial character.
            // All characters of the character array are accessed by
            // incrementing *lpt.
            
            for(int j=0; j<5; j++){
                IsoLabel[j] = *(lpt+j);               // Class field
                isoLabel[isoindex][j] = IsoLabel[j];  // Array
            }
            
        }
        
        // Partition function for this species (24 entries)
        
        void setpf(int i, double pfvalue){pf[i] = pfvalue; };
        
        void setfminus(double f, double FminusSum[]){
            fminus = f;                 // Class field
            FminusSum[isoindex] = f;    // Array
        }
        
        void setfplus(double f, double FplusSum[]){
            fplus = f;                  // Class field
            FplusSum[isoindex] = f;     // Array entry
        }
        
        void setkeff(double f){ keff=f; }
        
        void setdYdt(double d, double dYDt[]){
            dYdt = d; 
            dXdt = d*(double)A;
            dYDt[isoindex] = d;
        }
        
        void setdXdt(double d){
            dXdt = d;
            dYdt = d/(double)A;
        }
        
        
        // Public getter functions to return values of class fields
        
        int getZ() {return ZZ; };
        int getN() {return NN; };
        int getA() {return A; };
        double getY0() { return YY0; }
        double getY() {return YY; };
        double getX() {return XX; };
        double getM(){return MassExcess; };
        
        /* Overloaded getLabel function: getLabel() returns pointer to isoLabel character
         array; getLabel(k) returns kth character of isoLabel character array.
         
         Example usage:
        
            Species object1;
              .
              .
              .
            printf("\nLabel=%s", object1.getLabel());
          
            char label[5];
            for(k=0; k<5; k++){
               label[k] = object1.getLabel(k);
            }
            printf("\nLabel=%s", label);
            
         assuming that the data fields of object1 have already been populated using
         the setX functions.
        */
        
        char *getLabel() {return IsoLabel; };         // return pointer to label array
        char getLabel(int k) {return IsoLabel[k]; };  // return kth character of array
        
        // Partition function table
        
        double getpf(int i){return pf[i]; }
        
        // Partition function table temperatures 
        
        double getTpf(int i){return Tpf[i]; }
        
        double getfminus(){return fminus; }
        
        double getfplus(){return fplus; }
        
        double getkeff(){ return keff; }
        
        double getdYdt(){return dYdt; }
        
        double getdXdt(){return dXdt; }
        

    // Make data fields private, with external access to them through public setter 
    // and getter functions
    
    private:
        
        int isoindex;        // Index for species in isotope vector
        char IsoLabel[5];    // symbol for isotope
        int ZZ;              // proton number
        int NN;              // neutron number
        int A;               // atomic mass number
        double YY0;          // Abundance Y at beginning of timestep
        double YY;           // current abundance  Y = X/A
        double XX;           // mass fraction  X = Y*A
        double MassExcess;   // mass excess
        std::array<double, PF> pf;       // partition function entries
        double fplus;        // Total flux currently adding to abundance of this isotope
        double fminus;       // Total flux currently adding to abundance of this isotope
        double keff;         // Effective decay constant = fminus/YY
        double dYdt;         // Current dY/dt for this isotope
        double dXdt;         // Current dX/dt for this isotope
        
        // Temperatures in units of 10^9 K for partition function table (see pf[]). 
        
        const double Tpf[PF] = { 0.1f, 0.15f, 0.2f, 0.3f, 0.4f, 0.5f, 
            0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 3.5f, 
            4.0f, 4.5f, 5.0f, 6.0f, 7.0f, 8.0f, 9.0f, 10.0f };
        
};      // End of class Species
