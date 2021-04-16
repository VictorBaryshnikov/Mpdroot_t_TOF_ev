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

//максимум из двух чисел типа Double_t
Double_t MAX_2(Double_t Pdg_Mass, Double_t Hypotez_Mass)
{
if (Pdg_Mass >= Hypotez_Mass) return Pdg_Mass;
else return Hypotez_Mass;
}

// функция сортировки двумерного массива для сметченных трэков в выбранном ивенте [i][0] = индекс сметченного трэка, [i][1] = p - импульс
void qsotr_Data_n(Double_t** Data_n, Int_t left, Int_t right)
{
    Int_t i = left, j = right;
    Double_t temp, pivot = Data_n[ (left+right)/2 ][1];

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
            }

            i++; j--;
        }

    };

    if (left < j) qsotr_Data_n(Data_n, left, j);
    if (i < right) qsotr_Data_n(Data_n, i, right);
}

//функция подсчета погрешности импульса в зависимости от его значения, информация взята из TpcTdr-v07.pdf
Double_t Sigma_p( Double_t p) 
{
   Double_t dp_per_p;
   if (p <= 0.3) dp_per_p = 6.82 - 65.057936*p + 251.761904*p*p - 324.444444*p*p*p;
   else
       {
       if(p <= 1.2) dp_per_p = 1.3211148*p + 0.6804307;
       else dp_per_p = 1.325*p + 0.74;
       }

return(dp_per_p*p/100); //возвращаем уже сигму импульса
}

//функция перевода одномерного массива гипотез масс n треков в число: [0][1][1][0][1] -> 1101  
Int_t int_hipotez (Int_t* hipotez, Int_t number_elem_interval)
{
   Int_t result = 0;
   for (Int_t k = (number_elem_interval - 1); k >= 0; k--)
   {
      result = result + hipotez[k] * pow(10, (number_elem_interval - 1 - k) );
   }
   return (result);
}

//функция приращения гипотезы масс n треков 
void increase_elem_arr (Int_t* hipotez, Int_t number_elem_interval) //[0][0][0][0] -> [0][0][0][1] -> [0][0][1][0] -> [0][0][1][1] -> [0][1][0][0] -> [0][1][0][1] -> [0][1][1][0] -> [0][1][1][1] -> [1][0][0][0] -> [1][0][0][1] -> [1][0][1][0] -> [1][0][1][1] -> [1][1][0][0]   -> [1][1][0][1] -> [1][1][1][0] -> [1][1][1][1]
{

   hipotez[number_elem_interval - 1] = hipotez[number_elem_interval - 1] + 1;
   for (Int_t k = (number_elem_interval - 1); k > 0; k--)
   {
      if (hipotez[k] > 1)
      {
         hipotez[k-1] = hipotez[k-1] + 1;
         hipotez[k] = 0;
      }
   }
}

//функция вывода гипотезы масс n треков
void printf_hipotez (Int_t* hipotez, Int_t number_elem_interval)
{
   printf("hipotez: ");
   for (Int_t k = 0; k < number_elem_interval; k++)
   {
      printf("[%d] ", hipotez[k]);
   }
   printf(";\n ");
}

//
void printf_result (Double_t** result_n, Int_t N_intervals)
{
   //result_n[j][0] = time_tof_ev; // t_best_Ev
   //result_n[j][1] = time_tof_ev_sigma; // t_best_Ev_sigma
   //result_n[j][2] = Xi2; // Xi2
   //result_n[j][3] = doub_hipotez(hipotez, number_elem_interval); // гипотеза масс в виде числа '1101'
   //result_n[j][4] = num

   printf("\nresult: \n");
   for (Int_t j = 0; j < N_intervals; j++)
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

void record(int m[], int k[], int n) 
{
   for (int i = 0; i < n; i++)
   {
      k[i] = m[i];
   }
}

//заглушки для нейтарльных частиц, для ионов, для фотонов
Int_t plug_for_efficiency (Int_t Pdg_Code, int* plug)
{
   Int_t res = 1;
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
Int_t plug_motherid (Int_t MotherId, int plug)
{//Int_t MotherId = Int_t(Data_n[n_tracks][4]);

   Int_t res = 1;
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

void hipotez_null(Int_t number_elem_interval, Int_t* hipotez)
{
   for (Int_t k = 0; k < number_elem_interval; k++)
   {
      hipotez[k] = 0;
   }
}

//......................................................        class TSIGMA
class TSIGMA
{
private:

public:

	Int_t N_tofmatching; // N_tofmatching кол-во сметченных хитов в ивенте
        Double_t time_tof_ev_sigma; // t_best_Ev_sigma, ps

TSIGMA() // дефолтный конструктор 
	{	
                N_tofmatching = 0;	
		time_tof_ev_sigma = 0;
	}

TSIGMA(Int_t a0, Double_t a1) // конструктор от 2 чисел
	{	
		N_tofmatching = a0;
		time_tof_ev_sigma = a1;
	}

TSIGMA (const TSIGMA& tr) // Constructor for result
	{		
	        N_tofmatching = tr.N_tofmatching;
		time_tof_ev_sigma = tr.time_tof_ev_sigma;
	}

~TSIGMA()
	{
		//cout << "TSIGMA earsed." << endl;
	}

//......................................................      = 
TSIGMA operator=(const TSIGMA& a) 
	{
                if ( this == &a ) return *this; // проверка что это не S = S;
	        N_tofmatching = a.N_tofmatching;
		time_tof_ev_sigma = a.time_tof_ev_sigma;
                return *this;
	}
//......................................................      <<
friend ostream& operator<<(ostream& os, const TSIGMA& X)
	{
		os << "\nN_tofmatching = " << X.N_tofmatching << "; " << "\ntime_tof_ev_sigma = " << X.time_tof_ev_sigma << ", ps; " << endl;
		return os;
	}

//........................................... swap
void swap (TSIGMA& Y)
        {
                TSIGMA temp(*this);
                (*this) = Y;
                Y = temp;
        }

//........................................... set
void set (Int_t a0, Double_t a1)
        {
                N_tofmatching = a0;
		time_tof_ev_sigma = a1;
        }

//......................................................      +
TSIGMA operator +(const TSIGMA& Y)
	{
	        N_tofmatching = N_tofmatching + Y.N_tofmatching;
		time_tof_ev_sigma = time_tof_ev_sigma + Y.time_tof_ev_sigma;
                return (*this);
        }

//......................................................      -
TSIGMA operator /(const Int_t& X)
	{
                N_tofmatching = N_tofmatching / X;
		time_tof_ev_sigma = time_tof_ev_sigma / X;
                return (*this);
        }

};

// функция сортировки одномерного массива TSIGMA_n
void qsotr_TSIGMA_n(vector<TSIGMA> TSIGMA_n, Int_t left, Int_t right)
{
    Int_t i = left, j = right;
    Double_t pivot = TSIGMA_n[ int((left+right)/2) ].N_tofmatching;

    while (i <= j)
    {
        while (TSIGMA_n[i].N_tofmatching < pivot) i++;
        while (TSIGMA_n[j].N_tofmatching > pivot) j--;

        if (i <= j)
        {
            if (TSIGMA_n[i].N_tofmatching > TSIGMA_n[j].N_tofmatching)
            {
                TSIGMA_n[i].swap(TSIGMA_n[j]);
            }

            i++; j--;
        }

    };

    if (left < j) qsotr_TSIGMA_n(TSIGMA_n, left, j);
    if (i < right) qsotr_TSIGMA_n(TSIGMA_n, i, right);
}


//......................................................         class Result
class Result
{
private:
/*
	Double_t time_tof_ev; // t_best_Ev
        Double_t time_tof_ev_sigma; // t_best_Ev_sigm
        Double_t Xi2; // Xi2
        Int_t int_hipotez; // гипотеза масс в виде числа '1101'
        Int_t number_elem_interval; //кол.во элементов в интервале
*/  
public:

	Double_t time_tof_ev; // t_best_Ev
        Double_t time_tof_ev_sigma; // t_best_Ev_sigm
        Double_t Xi2; // Xi2
        Int_t int_hipotez; // гипотеза масс в виде числа '1101'
        Int_t number_elem_interval; //кол.во элементов в интервале

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

Result(Double_t a0, Double_t a1, Double_t a2, Int_t a3,  Int_t a4) // конструктор от 5 чисел 
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
void printf_xi2(Int_t Index_interval, Double_t efficiency_hypotez_Xi2min)
        {
                printf("Index interval %d; Xi2min = %f; efficiency_hypotez_Xi2min %f, '%'; \n", Index_interval, (*this).Xi2, efficiency_hypotez_Xi2min );
        }

//функция перевода числа int (100101) и кол.во элементов N в массив hypotez[N]
void int_in_hipotez( Int_t* hipot ) //1) гипотеза в виде 1101, 2) кол.во элементов в интервале 3) одномерный массив интов для гипотезы масс в виде 0 и 1
{//начало 0
   
   Int_t inumber_elem_interval = (*this).number_elem_interval;


   Int_t digit_number = 0; //разряд int_hipotez
   Int_t num = (*this).int_hipotez;

   if(num == 0 || num == 1)
   {
      digit_number = 1;
   }
   else
   {
      for (Int_t k = 2; k < 14; k++)
      {
         if ( num/pow(10, (k-1) ) >= 1 )
         {
            digit_number = k;
         }
      }
   }
   
   //printf("number = %d; digit_number = %d; inumber_elem_interval = %d;\n ", num, digit_number, inumber_elem_interval);
   Int_t n_lost_numbers = inumber_elem_interval - digit_number; //кол.во "потерянных" 0 при записей массива hipotez в число 
/*
   vector<int> vec;

   if (n_lost_numbers > 0)
   {
      for (int i = 0; i < n_lost_numbers; i++)
      {
         vec.push_back(0);
      }
   }

   if(digit_number == 1)
   {
      vec.push_back(num);
   }
   else
   {
      int del = digit_number - 1;
      for (int i = n_lost_numbers; i < inumber_elem_interval; i++)
      {
         vec.push_back( int(num/pow(10, del)) );

         if (vec[i] != 0 && del > 0 )
         {
            num = num - pow(10, del);
         }
         if ( del > 0 )
         {
            del = del - 1;
         }
      }
   }

   record(vec.data(), hipot, inumber_elem_interval);
*/

   if (n_lost_numbers > 0)
   {
      for (Int_t i = 0; i < n_lost_numbers; i++)
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
      Int_t del = digit_number - 1;
      for (Int_t i = n_lost_numbers; i < inumber_elem_interval; i++)
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

void printf_result (Result* result_n, Int_t N_intervals)
	{
	        printf("\nresult: \n");
                for (Int_t j = 0; j < N_intervals; j++)
                {
                   printf("[interval_index: %d]\n", j);
                   cout << result_n[j] << endl;
                }
                printf("\n ");
	}
//------------------------------------------------------------------------------------------------------------------------
#endif
