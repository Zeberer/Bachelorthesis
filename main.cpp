#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "array"
#include <algorithm>



using namespace std;

//Arrays zum speichern der einzelnen Variablen festlegen, sind immer die ersten Ableitungen gemeint

array<double, 100000> t1_array;
array<double, 100000> t2_array;
array<double, 100000> r1_array;
array<double, 100000> r2_array;
array<double, 100000> theta1_array;
array<double, 100000> theta2_array;
array<double, 100000> phi1_array;
array<double, 100000> phi2_array;

double M = 1;                   //Masse Schwarzes Loch auf eine Sonnenmasse festgelegt (Masse Testteilchen als null angenommen)
double G = 1;                   //Gravitationskonstante festgesetzt
double c = 1;                   //Lichtgeschwindigkeit festgesetzt


int sgn(double j) {
    if (j>0)
        return 1;
    else if (j<0)
        return -1;
    else
        return 0;
}





//8 Gleichungen als Funktionen definieren

double f1(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return t2;
}

double f2(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return r2;
}

double f3(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return theta2;
}

double f4(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return phi2;
}


double f5(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return - (M*G)/(pow(r1,2)*pow(c,2) - 2*M*G*r1) *t2*r2;
}

double f6(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return - (M*G*(1-((2*M*G)/(r1*pow(c,2)))))/(pow(r1,2)) * pow(c,2)* pow(t2,2)
           + ((M*G)/(pow(r1,2)*pow(c,2) - 2*M*G*r1)) * pow(r2,2)
           + (r1-((2*M*G)/(pow(c,2)))) * (pow(theta2,2) + pow(sin(theta1),2)*pow(phi2,2));
}

double f7(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return -(1/r1)*r2*theta2 + cos(theta1)*sin(theta1)* pow(phi2,2);
}

double f8(double t1, double t2, double r1, double r2, double theta1, double theta2, double phi1, double phi2, double tau) {
    return -(1/r1)*r2*phi2 - (cos(theta1)/sin(theta1))*theta2*phi2;
}









/*double f(double y, double t){          //f(double t, double r, double theta, double phi)
    return cos(y);                                 //hier alle geodätengleichungen angeben mit double y8 und dann gleichungen für alle y angeben
} */                                             //bzw geodätengleichung allg angeben und dann ausrechnen lassen
                                                //extra programm für berechnung christoffelsymbole
                                                //-> brauch ich nicht wenn ich geodätengleichung später mit euler-lagrange berechne


int main() {

    //Anfangswerte in kartesischen Koordinaten angeben und dann in Kugelkoordinaten umrechnen

    double x = 3;
    double y = 0;
    double z = 0;
    double dx = 0;
    double dy = 0;
    double dz = 0;

    //Orte

    t1_array[0] = 0;
    r1_array[0] = sqrt(pow(x, 2)+ pow(y,2) + pow(z, 2));
    cout << "r1 = " << r1_array[0] << endl;
    theta1_array[0] = acos(z/ r1_array[0]);
    if (x>0) {
       phi1_array[0] = atan(y/x);
    }

    if (x==0) {
        phi1_array[0] = sgn(y) * M_PI_2;
    }

    if (x<0 && y>=0) {
        phi1_array[0] = atan(y/x) + M_PI;
    }

    if (x<0 && y<0) {
        phi1_array[0] = atan(y/x) - M_PI;
    }

    //Geschwindigkeiten


    /*
    cout << "x = " << x << endl;

    double Z = x*dx+y*dy+z*dz;
    cout << "Z = " << Z << endl;
    double N = pow(x, 2) + pow(y,2) + pow(z, 2);
    cout << "N = " << N << endl;
    N = sqrt(N);
    cout << "N = " << N << endl;
     */

    r2_array[0] = (x*dx+y*dy+z*dz)/ sqrt(pow(x, 2) + pow(y,2) + pow(z, 2));
    cout << "r2 =" << r2_array[0] << endl;
    theta2_array[0] = sqrt(pow(x, 2)+ pow(y, 2)+ pow(z, 2)) * (pow(y, 2)+ pow(x, 2))*dz - ((x*dx+y*dy)*z)/(sqrt(pow(x, 2)+ pow(y, 2)) * (pow(x, 2)+ pow(y, 2)+ pow(z, 2)));
    if (x>0) {
        phi2_array[0] = sqrt(pow(x, 2)+ pow(y, 2)+ pow(z, 2))* sin(theta1_array[0]) * (x*dy-dx*y)/(pow(x, 2)+ pow(y, 2));
    }

    if (x==0) {
        phi2_array[0] = 0;
    }

    if (x<0 && y>=0) {
        phi2_array[0] = sqrt(pow(x, 2)+ pow(y, 2)+ pow(z, 2))* sin(theta1_array[0]) *(x*dy-dx*y)/(pow(x, 2)+ pow(y, 2));
    }

    if (x<0 && y<0) {
        phi2_array[0] = sqrt(pow(x, 2)+ pow(y, 2)+ pow(z, 2))* sin(theta1_array[0]) *(x*dy-dx*y)/(pow(x, 2)+ pow(y, 2));
    }




    //t1_array[0] = 0;
    //r1_array[0] = 100;
    //r2_array[0] = -0.1;
    //theta1_array[0] = 90 * M_PI/180;
    //theta2_array[0] = 0;
    //phi1_array[0] = 70 * M_PI/180;
    //phi2_array[0] = 0.05;


    double A = pow(c, 2)*(1-((2*M*G)/(pow(c,2)*r1_array[0])));
    double B = 1/(1-((2*M*G)/(r1_array[0]* pow(c, 2))));

    cout<< "A =" << A << endl;
    cout<< "B =" << B << endl;


    double t2 = -((-B* pow(r2_array[0], 2) - pow(r1_array[0], 2)* pow(theta2_array[0], 2) - pow(r1_array[0], 2)* pow(sin(theta1_array[0]), 2)* pow(phi2_array[0], 2))/A);
    //cout << "t2 = " << t2 << endl;
    t2 = sqrt(t2);
    t2_array[0] = t2;
    cout << "t2 = " << t2 << endl;





    //Gleichung für die Normierung, muss immer null sein
    //double N = A* pow(t2_array[0], 2) - B* pow(r2_array[0], 2) - pow(r1_array[0], 2)* pow(theta2_array[0], 2) - pow(r1_array[0], 2)*pow(sin(theta1_array[0]), 2)*pow(phi2_array[0], 2);
    //cout << "N = " << N << endl;





    double tau = 0;             //Startwert für tau festlegen (richtige Zeit/Eigenzeit)
    double dtau = 0.1;         //Zeitschritt den es immer weiter geht festlegen
    double tauend = 5000;        //Endwert für tau festlegen, liegt der Wert über 2.46 dann gibt es nen segmentation fault,
                                // heißt der Speicher ist zu voll (array größer gemacht damit es geht)

    //k variablen deklarieren und auf gemeinsamen startwert initialisieren

    double k1f1 = 0, k2f1 = 0, k3f1 = 0, k4f1 = 0;
    double k1f2 = 0, k2f2 = 0, k3f2 = 0, k4f2 = 0;
    double k1f3 = 0, k2f3 = 0, k3f3 = 0, k4f3 = 0;
    double k1f4 = 0, k2f4 = 0, k3f4 = 0, k4f4 = 0;
    double k1f5 = 0, k2f5 = 0, k3f5 = 0, k4f5 = 0;
    double k1f6 = 0, k2f6 = 0, k3f6 = 0, k4f6 = 0;
    double k1f7 = 0, k2f7 = 0, k3f7 = 0, k4f7 = 0;
    double k1f8 = 0, k2f8 = 0, k3f8 = 0, k4f8 = 0;

    double t2_berechnet;

    FILE *eingabedatei;
    eingabedatei = fopen("C:\\Users\\SoZeb\\PycharmProjects\\pythonProject\\Gerade.csv", "w");


    //Variable für den Punkt im Array initialisieren

    int i = 0;


    //Werte vor Beginn der Schleife ausgeben (Startwerte)

    printf("%e %e %e %e %e %e %e %e %e\n", t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
    //fprintf(eingabedatei, "%e %e %e %e %e %e %e %e\n", X1, Y1, Z1, X2, Y2, Z2, dtau, tauend);
    fprintf(eingabedatei, "%e %e %e\n",  r1_array[i], theta1_array[i], phi1_array[i]);




    //Schleife die von tau =1 bis zur Größe des Arrays immer einen Zeitschritt weitergeht

    while(tau<tauend){


        //Koeffizienten für die erste Gleichung berechnen und zu neuem Punkt aufaddieren

        k1f1 = f1(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f1 = f1(t1_array[i]+dtau*(k1f1/2.), t2_array[i]+dtau*(k1f1/2.), r1_array[i]+dtau*(k1f1/2.), r2_array[i]+dtau*(k1f1/2.), theta1_array[i]+dtau*(k1f1/2.), theta2_array[i]+dtau*(k1f1/2.), phi1_array[i]+dtau*(k1f1/2.), phi2_array[i]+dtau*(k1f1/2.), tau + dtau/2.);
        k3f1 = f1(t1_array[i]+dtau*(k2f1/2.), t2_array[i]+dtau*(k2f1/2.), r1_array[i]+dtau*(k2f1/2.), r2_array[i]+dtau*(k2f1/2.), theta1_array[i]+dtau*(k2f1/2.), theta2_array[i]+dtau*(k2f1/2.), phi1_array[i]+dtau*(k2f1/2.), phi2_array[i]+dtau*(k2f1/2.), tau + dtau/2.);
        k4f1 = f1(t1_array[i]+dtau*k3f1, t2_array[i]+dtau*k3f1, r1_array[i]+dtau*k3f1, r2_array[i]+dtau*k3f1, theta1_array[i]+dtau*k3f1, theta2_array[i]+dtau*k3f1, phi1_array[i]+dtau*k3f1, phi2_array[i]+dtau*k3f1, tau + dtau);

        t1_array[i+1] = t1_array[i] + (dtau/6.) * (k1f1 + 2.*k2f1 + 2.*k3f1 + k4f1);
        //printf("t1: %e\n", t1_array[i+1.]);



        k1f2 = f2(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f2 = f2(t1_array[i]+dtau*(k1f2/2.), t2_array[i]+dtau*(k1f2/2.), r1_array[i]+dtau*(k1f2/2.), r2_array[i]+dtau*(k1f2/2.), theta1_array[i]+dtau*(k1f2/2.), theta2_array[i]+dtau*(k1f2/2.), phi1_array[i]+dtau*(k1f2/2.), phi2_array[i]+dtau*(k1f2/2.), tau + dtau/2.);
        k3f2 = f2(t1_array[i]+dtau*(k2f2/2.), t2_array[i]+dtau*(k2f2/2.), r1_array[i]+dtau*(k2f2/2.), r2_array[i]+dtau*(k2f2/2.), theta1_array[i]+dtau*(k2f2/2.), theta2_array[i]+dtau*(k2f2/2.), phi1_array[i]+dtau*(k2f2/2.), phi2_array[i]+dtau*(k2f2/2.), tau + dtau/2.);
        k4f2 = f2(t1_array[i]+dtau*k3f2, t2_array[i]+dtau*k3f2, r1_array[i]+dtau*k3f2, r2_array[i]+dtau*k3f2, theta1_array[i]+dtau*k3f2, theta2_array[i]+dtau*k3f2, phi1_array[i]+dtau*k3f2, phi2_array[i]+dtau*k3f2, tau + dtau);

        r1_array[i+1] = r1_array[i] + (dtau/6.) * (k1f2 + 2.*k2f2 + 2.*k3f2 + k4f2);
        //printf("r1: %e\n", r1_array[i+1.]);



        k1f3 = f3(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f3 = f3(t1_array[i]+dtau*(k1f3/2.), t2_array[i]+dtau*(k1f3/2.), r1_array[i]+dtau*(k1f3/2.), r2_array[i]+dtau*(k1f3/2.), theta1_array[i]+dtau*(k1f3/2.), theta2_array[i]+dtau*(k1f3/2.), phi1_array[i]+dtau*(k1f3/2.), phi2_array[i]+dtau*(k1f3/2.), tau + dtau/2.);
        k3f3 = f3(t1_array[i]+dtau*(k2f3/2.), t2_array[i]+dtau*(k2f3/2.), r1_array[i]+dtau*(k2f3/2.), r2_array[i]+dtau*(k2f3/2.), theta1_array[i]+dtau*(k2f3/2.), theta2_array[i]+dtau*(k2f3/2.), phi1_array[i]+dtau*(k2f3/2.), phi2_array[i]+dtau*(k2f3/2.), tau + dtau/2.);
        k4f3 = f3(t1_array[i]+dtau*k3f3, t2_array[i]+dtau*k3f3, r1_array[i]+dtau*k3f3, r2_array[i]+dtau*k3f3, theta1_array[i]+dtau*k3f3, theta2_array[i]+dtau*k3f3, phi1_array[i]+dtau*k3f3, phi2_array[i]+dtau*k3f3, tau + dtau);

        theta1_array[i+1] = theta1_array[i] + (dtau/6.) * (k1f3 + 2.*k2f3 + 2.*k3f3 + k4f3);
        //printf("theta1: %e\n", theta1_array[i+1.]);



        k1f4 = f4(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f4 = f4(t1_array[i]+dtau*(k1f4/2.), t2_array[i]+dtau*(k1f4/2.), r1_array[i]+dtau*(k1f4/2.), r2_array[i]+dtau*(k1f4/2.), theta1_array[i]+dtau*(k1f4/2.), theta2_array[i]+dtau*(k1f4/2.), phi1_array[i]+dtau*(k1f4/2.), phi2_array[i]+dtau*(k1f4/2.), tau + dtau/2.);
        k3f4 = f4(t1_array[i]+dtau*(k2f4/2.), t2_array[i]+dtau*(k2f4/2.), r1_array[i]+dtau*(k2f4/2.), r2_array[i]+dtau*(k2f4/2.), theta1_array[i]+dtau*(k2f4/2.), theta2_array[i]+dtau*(k2f4/2.), phi1_array[i]+dtau*(k2f4/2.), phi2_array[i]+dtau*(k2f4/2.), tau + dtau/2.);
        k4f4 = f4(t1_array[i]+dtau*k3f4, t2_array[i]+dtau*k3f4, r1_array[i]+dtau*k3f4, r2_array[i]+dtau*k3f4, theta1_array[i]+dtau*k3f4, theta2_array[i]+dtau*k3f4, phi1_array[i]+dtau*k3f4, phi2_array[i]+dtau*k3f4, tau + dtau);

        phi1_array[i+1] = phi1_array[i] + (dtau/6.) * (k1f4 + 2.*k2f4 + 2.*k3f4 + k4f4);
        //printf("phi1: %e\n", phi1_array[i+1.]);



        k1f5 = f5(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f5 = f5(t1_array[i]+dtau*(k1f5/2.), t2_array[i]+dtau*(k1f5/2.), r1_array[i]+dtau*(k1f5/2.), r2_array[i]+dtau*(k1f5/2.), theta1_array[i]+dtau*(k1f5/2.), theta2_array[i]+dtau*(k1f5/2.), phi1_array[i]+dtau*(k1f5/2.), phi2_array[i]+dtau*(k1f5/2.), tau + dtau/2.);
        k3f5 = f5(t1_array[i]+dtau*(k2f5/2.), t2_array[i]+dtau*(k2f5/2.), r1_array[i]+dtau*(k2f5/2.), r2_array[i]+dtau*(k2f5/2.), theta1_array[i]+dtau*(k2f5/2.), theta2_array[i]+dtau*(k2f5/2.), phi1_array[i]+dtau*(k2f5/2.), phi2_array[i]+dtau*(k2f5/2.), tau + dtau/2.);
        k4f5 = f5(t1_array[i]+dtau*k3f5, t2_array[i]+dtau*k3f5, r1_array[i]+dtau*k3f5, r2_array[i]+dtau*k3f5, theta1_array[i]+dtau*k3f5, theta2_array[i]+dtau*k3f5, phi1_array[i]+dtau*k3f5, phi2_array[i]+dtau*k3f5, tau + dtau);

        t2_array[i+1] = t2_array[i] + (dtau/6.) * (k1f5 + 2.*k2f5 + 2.*k3f5 + k4f5);
        //printf("t2: %e\n", t2_array[i+1.]);



        k1f6 = f6(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f6 = f6(t1_array[i]+dtau*(k1f6/2.), t2_array[i]+dtau*(k1f6/2.), r1_array[i]+dtau*(k1f6/2.), r2_array[i]+dtau*(k1f6/2.), theta1_array[i]+dtau*(k1f6/2.), theta2_array[i]+dtau*(k1f6/2.), phi1_array[i]+dtau*(k1f6/2.), phi2_array[i]+dtau*(k1f6/2.), tau + dtau/2.);
        k3f6 = f6(t1_array[i]+dtau*(k2f6/2.), t2_array[i]+dtau*(k2f6/2.), r1_array[i]+dtau*(k2f6/2.), r2_array[i]+dtau*(k2f6/2.), theta1_array[i]+dtau*(k2f6/2.), theta2_array[i]+dtau*(k2f6/2.), phi1_array[i]+dtau*(k2f6/2.), phi2_array[i]+dtau*(k2f6/2.), tau + dtau/2.);
        k4f6 = f6(t1_array[i]+dtau*k3f6, t2_array[i]+dtau*k3f6, r1_array[i]+dtau*k3f6, r2_array[i]+dtau*k3f6, theta1_array[i]+dtau*k3f6, theta2_array[i]+dtau*k3f6, phi1_array[i]+dtau*k3f6, phi2_array[i]+dtau*k3f6, tau + dtau);

        r2_array[i+1] = r2_array[i] + (dtau/6.) * (k1f6 + 2.*k2f6 + 2.*k3f6 + k4f6);
        //printf("r2: %e\n", r2_array[i+1.]);



        k1f7 = f7(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        k2f7 = f7(t1_array[i]+dtau*(k1f7/2.), t2_array[i]+dtau*(k1f7/2.), r1_array[i]+dtau*(k1f7/2.), r2_array[i]+dtau*(k1f7/2.), theta1_array[i]+dtau*(k1f7/2.), theta2_array[i]+dtau*(k1f7/2.), phi1_array[i]+dtau*(k1f7/2.), phi2_array[i]+dtau*(k1f7/2.), tau + dtau/2.);
        k3f7 = f7(t1_array[i]+dtau*(k2f7/2.), t2_array[i]+dtau*(k2f7/2.), r1_array[i]+dtau*(k2f7/2.), r2_array[i]+dtau*(k2f7/2.), theta1_array[i]+dtau*(k2f7/2.), theta2_array[i]+dtau*(k2f7/2.), phi1_array[i]+dtau*(k2f7/2.), phi2_array[i]+dtau*(k2f7/2.), tau + dtau/2.);
        k4f7 = f7(t1_array[i]+dtau*k3f7, t2_array[i]+dtau*k3f7, r1_array[i]+dtau*k3f7, r2_array[i]+dtau*k3f7, theta1_array[i]+dtau*k3f7, theta2_array[i]+dtau*k3f7, phi1_array[i]+dtau*k3f7, phi2_array[i]+dtau*k3f7, tau + dtau);

        theta2_array[i+1] = theta2_array[i] + (dtau/6.) * (k1f7 + 2.*k2f7 + 2.*k3f7 + k4f7);
        //printf("theta2: %e\n", theta2_array[i+1.]);



        k1f8 = f8(t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);     //Wert ist zu klein, gibt halbiert wahrscheinlich einfach null
        k2f8 = f8(t1_array[i]+dtau*(k1f8/2.), t2_array[i]+dtau*(k1f8/2.), r1_array[i]+dtau*(k1f8/2.), r2_array[i]+dtau*(k1f8/2.), theta1_array[i]+dtau*(k1f8/2.), theta2_array[i]+dtau*(k1f8/2.), phi1_array[i]+dtau*(k1f8/2.), phi2_array[i]+dtau*(k1f8/2.), tau + dtau/2.);
        k3f8 = f8(t1_array[i]+dtau*(k2f8/2.), t2_array[i]+dtau*(k2f8/2.), r1_array[i]+dtau*(k2f8/2.), r2_array[i]+dtau*(k2f8/2.), theta1_array[i]+dtau*(k2f8/2.), theta2_array[i]+dtau*(k2f8/2.), phi1_array[i]+dtau*(k2f8/2.), phi2_array[i]+dtau*(k2f8/2.), tau + dtau/2.);
        k4f8 = f8(t1_array[i]+dtau*k3f8, t2_array[i]+dtau*k3f8, r1_array[i]+dtau*k3f8, r2_array[i]+dtau*k3f8, theta1_array[i]+dtau*k3f8, theta2_array[i]+dtau*k3f8, phi1_array[i]+dtau*k3f8, phi2_array[i]+dtau*k3f8, tau + dtau);

        phi2_array[i+1] = phi2_array[i] + (dtau/6.) * (k1f8 + 2.*k2f8 + 2.*k3f8 + k4f8);
        //printf("phi2: %e\n", phi2_array[i+1.]);


       if(r1_array[i] < 2){
           break;
           cout << "Im Schwarzen Loch" << endl;
       }


      //double Geschwindigkeit = A * pow(c, 2) * pow(t2_array[i], 2) - B * pow(r2_array[i], 2) -
      //        pow(r1_array[i], 2) * pow(theta2_array[i], 2) - pow(r1_array[i], 2) * pow(sin(theta1_array[i]), 2) * pow(phi2_array[i], 2);

      //cout << "Geschwindigkeit = " << Geschwindigkeit << endl;





        //Variable zum festlegen des Punkts im Array eins erhöhen
        i ++;
        tau = tau + dtau;


        //Werte nach Berechnung des neuen Werts ausgeben
        printf("%e %e %e %e %e %e %e %e %e\n", t1_array[i], t2_array[i], r1_array[i], r2_array[i], theta1_array[i], theta2_array[i], phi1_array[i], phi2_array[i], tau);
        fprintf(eingabedatei, "%e %e %e \n", r1_array[i], theta1_array[i], phi1_array[i]);






    }





    fclose(eingabedatei);



    return 0;

}




