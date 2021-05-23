//------------------------------------------------------------------------------------------------------------------------
#ifndef __LIBRARY_ESTIMATION_H
#define __LIBRARY_ESTIMATION_H 1

#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;
//------------------------------------------------------------------------------------------------------------------------

//максимум из двух чисел
double MAX_2(double Pdg_Mass, double Hypotez_Mass)
{
if (Pdg_Mass >= Hypotez_Mass) return Pdg_Mass;
else return Hypotez_Mass;
}

//min из двух чисел
double MIN_2(double A, double B)
{
if (A >= B) return B;
else return A;
}

// функция сортировки двумерного массива для сметченных трэков в выбранном ивенте [i][0] = индекс сметченного трэка, [i][1] = p - импульс
void qsotr_Data_n(double** Data_n, int left, int right)
{
    int i = left, j = right;
    double temp, pivot = Data_n[ (left+right)/2 ][1];

    while (i <= j)
    {
        while (Data_n[i][1] < pivot) i++;
        while (Data_n[j][1] > pivot) j--;

        if (i <= j)
        {
            if (Data_n[i][1] > Data_n[j][1])
            {
                temp = Data_n[i][1]; Data_n[i][1] = Data_n[j][1]; Data_n[j][1] = temp;
                temp = Data_n[i][0]; Data_n[i][0] = Data_n[j][0]; Data_n[j][0] = temp;
                temp = Data_n[i][2]; Data_n[i][2] = Data_n[j][2]; Data_n[j][2] = temp;
                temp = Data_n[i][3]; Data_n[i][3] = Data_n[j][3]; Data_n[j][3] = temp;
                temp = Data_n[i][4]; Data_n[i][4] = Data_n[j][4]; Data_n[j][4] = temp;
                temp = Data_n[i][5]; Data_n[i][5] = Data_n[j][5]; Data_n[j][5] = temp;
                temp = Data_n[i][6]; Data_n[i][6] = Data_n[j][6]; Data_n[j][6] = temp;
            }

            i++; j--;
        }

    };

    if (left < j) qsotr_Data_n(Data_n, left, j);
    if (i < right) qsotr_Data_n(Data_n, i, right);
}

//функция подсчета погрешности импульса в зависимости от его значения, информация взята из TpcTdr-v07.pdf
double Sigma_p( double p) 
{
   double dp_per_p;
   if (p <= 0.3) dp_per_p = 6.82 - 65.057936*p + 251.761904*p*p - 324.444444*p*p*p;
   else
       {
       if(p <= 1.2) dp_per_p = 1.3211148*p + 0.6804307;
       else dp_per_p = 1.325*p + 0.74;
       }

return(dp_per_p*p/100); //возвращаем уже сигму импульса
}

//функция перевода одномерного массива гипотез масс n треков в число: [0][1][1][0][1] -> 1101  
long double int_hipotez (int* hipotez, int number_elem_interval)
{
   long double result = 0;
   for (int k = (number_elem_interval - 1); k >= 0; k--)
   {
      result = result + hipotez[k] * pow(10, (number_elem_interval - 1 - k) );
      //cout <<"\n" << result << endl;
   }
 
   return (result);
}

//функция приращения гипотезы масс n треков 
void increase_elem_arr (int* hipotez, int number_elem_interval) //[0][0][0][0] -> [0][0][0][1] -> [0][0][1][0] -> [0][0][1][1] -> [0][1][0][0] -> [0][1][0][1] -> [0][1][1][0] -> [0][1][1][1] -> [1][0][0][0] -> [1][0][0][1] -> [1][0][1][0] -> [1][0][1][1] -> [1][1][0][0]   -> [1][1][0][1] -> [1][1][1][0] -> [1][1][1][1]
{

   hipotez[number_elem_interval - 1] = hipotez[number_elem_interval - 1] + 1;
   for (int k = (number_elem_interval - 1); k > 0; k--)
   {
      if (hipotez[k] > 1)
      {
         hipotez[k-1] = hipotez[k-1] + 1;
         hipotez[k] = 0;
      }
   }
}

//функция вывода гипотезы масс n треков
void printf_hipotez (int* hipotez, int number_elem_interval)
{
   printf("hipotez: ");
   for (int k = 0; k < number_elem_interval; k++)
   {
      printf("[%d] ", hipotez[k]);
   }
   printf(";\n ");
}

//
void printf_result (double** result_n, int N_intervals)
{
   //result_n[j][0] = time_tof_ev; // t_best_Ev
   //result_n[j][1] = time_tof_ev_sigma; // t_best_Ev_sigma
   //result_n[j][2] = Xi2; // Xi2
   //result_n[j][3] = doub_hipotez(hipotez, number_elem_interval); // гипотеза масс в виде числа '1101'
   //result_n[j][4] = num

   printf("\nresult: \n");
   for (int j = 0; j < N_intervals; j++)
   {
      printf("[interval_index: %d][0]: t_best_Ev(time_tof_ev) = %f, ns;\n", j, result_n[j][0]);
      printf("[interval_index: %d][1]: t_best_Ev_sigma(time_tof_ev_sigma) = %f, ns;\n", j, result_n[j][1]);
      printf("[interval_index: %d][2]: Xi2 = %f;\n", j, result_n[j][2]);
      printf("[interval_index: %d][3]: hipotez_mass = %d;\n", j, Int_t(result_n[j][3]));
      printf("[interval_index: %d][4]: number_elem_interval = %d;\n", j, Int_t(result_n[j][4]));
      printf("\n ");
   }
   printf("\n ");
}

//заглушки для нейтарльных частиц, для ионов, для фотонов
int plug_for_efficiency (int Pdg_Code, int* plug)
{
   int res = 1;
   //int plug[3] = {0, 0, 0}; - включаем все заглушки
   if(plug[0] == 0) //заглушка для нейтральных частиц
   {
      //                pi0                    Xi0                       n                       D0                 Sigma0
      if( abs(Pdg_Code) == 111 || abs(Pdg_Code) == 311 || abs(Pdg_Code) == 2112 || abs(Pdg_Code) == 421 || abs(Pdg_Code) == 3212 ) 
      {
         //printf("plug for neutral particles; \n");
         res = 0;
         return res;
      }
      
   }

   if(plug[1] == 0) //заглушка для ионов
   {
      //          
      if( abs(Pdg_Code) > 10000) 
      {
         //printf("plug for ions; \n");
         res = 0;
         return res;
      }
      
   }

   if(plug[2] == 0) //заглушка для фотонов
   {
      //          
      if( abs(Pdg_Code) == 22) 
      {
         //printf("plug for photons; \n");
         res = 0;
         return res;
      }
      
   }
   
   return res;
}

//заглушка по MotherId
int plug_motherid (int MotherId, int plug)
{//int MotherId = Int_t(Data_n[n_tracks][4]);

   int res = 1;
   if(plug == 0)
   {
      if (MotherId != -1) 
      {
         res = 0;
         return res;
      }
   }
   return res;
}

void hipotez_null(int number_elem_interval, int* hipotez)
{
   for (int k = 0; k < number_elem_interval; k++)
   {
      hipotez[k] = 0;
   }
}

//......................................................         class Result
class Result
{
private:
/*
	double time_tof_ev; // t_best_Ev
        double time_tof_ev_sigma; // t_best_Ev_sigm
        double Xi2; // Xi2
        int int_hipotez; // гипотеза масс в виде числа '1101'
        int number_elem_interval; //кол.во элементов в интервале
*/  
public:

	double time_tof_ev; // t_best_Ev
        double time_tof_ev_sigma; // t_best_Ev_sigm
        double Xi2; // Xi2
        long double int_hipotez; // гипотеза масс в виде числа '1101'
        int number_elem_interval; //кол.во элементов в интервале

Result() // дефолтный конструктор 
	{	
                //cout << "\nDefault constructor" << endl;	
		time_tof_ev = 0;
		time_tof_ev_sigma = 0;
                Xi2 = 0;
                int_hipotez = 0; 
                number_elem_interval = 0;
		//cout << "result created." << endl;
	}

Result(double a0, double a1, double a2, long double a3,  int a4) // конструктор от 5 чисел 
	{
                //cout << "\nconstructor" << endl;	
		time_tof_ev = a0;
		time_tof_ev_sigma = a1;
                Xi2 = a2;
                int_hipotez = a3; 
                number_elem_interval = a4;
		//cout << "result created." << endl;
	}

Result (const Result& tr) // Constructor for result
	{
		//cout << "\nConstructor for Rational_number." << endl;
	        time_tof_ev = tr.time_tof_ev;
		time_tof_ev_sigma = tr.time_tof_ev_sigma;
		Xi2 = tr.Xi2;
		int_hipotez = tr.int_hipotez;
		number_elem_interval = tr.number_elem_interval;
		//cout << "result created." << endl;
	}

~Result()
	{
		//cout << "result earsed." << endl;
	}

//......................................................      = 
Result operator=(const Result& a) 
	{
                if ( this == &a ) return *this; // проверка что это не S = S;
	        time_tof_ev = a.time_tof_ev;
		time_tof_ev_sigma = a.time_tof_ev_sigma;
		Xi2 = a.Xi2;
		int_hipotez = a.int_hipotez;
                number_elem_interval = a.number_elem_interval;
        return *this;
	}
//......................................................      << (+)
friend ostream& operator<<(ostream& os, const Result& X)
	{
		os << "\ntime_tof_ev = " << X.time_tof_ev << ", ns; " << "\ntime_tof_ev_sigma = " << X.time_tof_ev_sigma << ", ns; " << "\nXi2 = " << X.Xi2 << "; " << "\nint_hipotez = " << X.int_hipotez << "; " << "\nnumber_elem_interval = " << X.number_elem_interval << "; " << endl;
		return os;
	}

//......................................................      printf_xi2 (+)
void printf_xi2(int Index_interval, double efficiency_hypotez_Xi2min)
        {
                printf("Index interval %d; Xi2min = %f; efficiency_hypotez_Xi2min %f, '%'; \n", Index_interval, (*this).Xi2, efficiency_hypotez_Xi2min );
        }

//функция перевода числа int (100101) и кол.во элементов N в массив hypotez[N]
void int_in_hipotez( int* hipot ) //1) гипотеза в виде 1101, 2) кол.во элементов в интервале 3) одномерный массив интов для гипотезы масс в виде 0 и 1
{//начало 0
   
   int inumber_elem_interval = (*this).number_elem_interval;


   int digit_number = 0; //разряд int_hipotez
   long double num = (*this).int_hipotez;

   //cout << "int_hipotez :" << int_hipotez << ", " << num << endl;

   if(num == 0 || num == 1)
   {
      digit_number = 1;
   }
   else
   {
      for (int k = 2; k < 22; k++)
      {
         if ( Int_t(num/pow(10, (k-1) ) ) >= 1 )
         {
            //cout << "Int_t(num/pow(10, (k-1) ) ) :" << Int_t(num/pow(10, (k-1) ) ) << endl;
            digit_number = k;
         }
      }
   }
   
   //printf("\n digit_number = %d; inumber_elem_interval = %d;\n ", digit_number, inumber_elem_interval);

   int n_lost_numbers = inumber_elem_interval - digit_number; //кол.во "потерянных" 0 при записей массива hipotez в число 

   if (n_lost_numbers > 0)
   {
      for (int i = 0; i < n_lost_numbers; i++)
      {
         hipot[i] = 0;
      }
   }
   
   if(digit_number == 1)
   {
      hipot[inumber_elem_interval - 1] = num;
   }
   else
   {
      int del = digit_number - 1;
      for (int i = n_lost_numbers; i < inumber_elem_interval; i++)
      {
         hipot[i] = Int_t(num/pow(10, del));

         if (hipot[i] != 0 && del > 0 )
         {
            num = num - pow(10, del);
         }
         if ( del > 0 )
         {
            del = del - 1;
         }
      }
   }

}//конец 0

};

void printf_result (Result* result_n, int N_intervals)
	{
	        printf("\nresult: \n");
                for (int j = 0; j < N_intervals; j++)
                {
                   printf("[interval_index: %d]\n", j);
                   cout << result_n[j] << endl;
                }
                printf("\n ");
	}


//......................................................         class AFP - Array for particle
class AFP //Array for particle
{
public:
        
	double t_exp;
	double sigma_pid;
        double p_t;

AFP() // дефолтный конструктор 
	{	
                //cout << "\nDefault constructor" << endl;	
		t_exp = 0;
                sigma_pid = 0;
                p_t = 0;
		//cout << "result created." << endl;
	}

AFP( double t, double sigma, double pt ) 
	{	
                //cout << "\n Constructor" << endl;	
		t_exp = t;
                sigma_pid = sigma;
                p_t = pt;
		//cout << "AFP created." << endl;
	}

//......................................................      << (+)
friend ostream& operator<<(ostream& os, const AFP& X)
	{
             //for (int i = 0; i < X.t_exp.size(); i++)
             //{
		os << "\n time_exp = " << X.t_exp << ", ns; " << "  Sigma_pid = " << X.sigma_pid << "; " << "   P_t = " << X.p_t << " GeV / c; " << endl;
             //}
		return os;
	}
};
//------------------------------------------------------------------------------------------------------------------------
#endif
