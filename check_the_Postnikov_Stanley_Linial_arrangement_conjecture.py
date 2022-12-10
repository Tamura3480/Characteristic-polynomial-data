# -*- coding: utf-8 -*-

#Python version 3.8.5
import math
from fractions import Fraction
from typing import List

from compute_the_characteristic_polynomial import (
    generate_generalized_Eulerian_polynomial,
    generate_Ehrhart_quasi_polynomial,
    generate_average_quasi_polynomial,
    shift_operator_act_on_polynomial,
    generate_c
)

#Compute the polynomial f(\sqrt{-1}t+a)
#input  :the polynomial f(t) and a
#output :the polynomial f(\sqrt{-1}t+a)
def shift_and_rot_polynomial(polynomial :List[int] ,a :int) -> List[int] :
    shift_operator=[0]*(a+1)
    shift_operator[a]=1
    d=len(polynomial)
    new_polynomial=shift_operator_act_on_polynomial(shift_operator,polynomial,a=-1)
    i=0
    b=3 if d%2==0 else 2
    while 4*i+b<d:
        new_polynomial[4*i+b]=-new_polynomial[4*i+b]
        i+=1
    return new_polynomial

#Compute the derivative of the polynomial f(t)
#input  :the polynomial f(t)
#output :the polynomial f'(t)
def differentiate_polynomial(polynomial :List[int]) -> List[int]:
    d=len(polynomial)
    new_polynomial=[0]*(d-1)
    for i in range(d-1):
        new_polynomial[i]=(i+1)*polynomial[i+1]
    return new_polynomial

#Compute the value of the polynomial at point a
#input  :the polynomial f(t)　and point a
#output :the value f(a)
def calc_value(polynomial :List[int] ,point :int) -> int:
    value=0
    d=len(polynomial)
    for i in range(d):
        value+=polynomial[i]*(point**(i))
    return value

#N(a):=the number of sign changes in the sequence {f(a),f'(a),f''(a), ...,f^((n))(a)}
#Compute N(start_point)-N(end_point) of f
#input  :the polynomial f(t) and the interval [start_point,end_point]
#output :N(start_point)-N(end_point)
def count_change_of_sign(polynomial :List[int] ,start_point :int ,end_point :int) -> int:
    d=len(polynomial)
    new_polynomial=polynomial[:]#list deep_copy
    start_value1=calc_value(new_polynomial,start_point)
    end_value1=calc_value(new_polynomial,end_point)
    if start_value1*end_value1!=0:
        start_count=0
        end_count=0
        for i in range(d-1):
            new_polynomial=differentiate_polynomial(new_polynomial)
            start_value2=calc_value(new_polynomial,start_point)
            end_value2=calc_value(new_polynomial,end_point)
            if start_value2!=0:
                if start_value1*start_value2<0:
                    start_count+=1
                start_value1=start_value2
            if end_value2!=0:
                if end_value1*end_value2<0:
                    end_count+=1
                end_value1=end_value2
    else:
        start_count=1 if start_value1==0 else 0
        end_count=0
    return start_count-end_count

#Compute the lower bound of the number of real roots of the polynomial f using the Fourier-Budan Theorem
#input  :the polynomial f and the interval [start_point, end_point]
#output :the lower bound of the number of real roots of polynomial f
def Fourier_Budan(polynomial :List[int] ,start_point :int ,end_point :int) -> int:
    lb_nsol=0
    d=len(polynomial)
    for i in range(start_point,end_point,1):
        Ndiff=count_change_of_sign(polynomial,i,i+1)
        lb_nsol=lb_nsol+1  if Ndiff==1 else lb_nsol
    return lb_nsol


#main program
#Check the Postnikov-Stanley Linial arrangement conjecture
#When we determine the variable "phi" and run this program, we can check the conjecture for the Linial arrangement of "phi".
#If the conjecture holds, then the sentence "All the roots of the characteristic polynomial have the same real part nh/2" appears.
#Only "E_6","E_7","E_8" or "F_4" can be selected as "phi"

#We transform the characteristic polynomial to the appropriate polynomial.
#Specifically, the characteristic polynomial is shifted by -nh/2 and then rotated 90 degrees around the origin.
#We count the number of real roots of the transformed characteristic polynomial using the Fourier Budan Theorem.
#If all the roots of the transformed characteristic polynomial are real, then all the roots of the characteristic polynomial have the same real part nh/2.

phi="E_8"

rad_period={"E_6":6,"E_7":6,"E_8":30,"F_4":6}
arrangement_list=list(set([math.gcd(i,rad_period[phi])-1 for i in range(rad_period[phi])]))#list of n of Linial arrangement \mathcal{A}^[1,n]
del arrangement_list[0]


generalized_Eulerian_polynomial=generate_generalized_Eulerian_polynomial(phi)#R_{phi}(t)
Ehrhart_quasi_polynomial=generate_Ehrhart_quasi_polynomial(phi)#L_{phi}(t)
period=len(Ehrhart_quasi_polynomial)
h=sum(generate_c(phi))+1

for n in arrangement_list:
    average_Ehrhart_quasi_polynomial=generate_average_quasi_polynomial(math.gcd(n+1,period),Ehrhart_quasi_polynomial)#\tilde{L}^{gcd(n+1,period)}_{phi}(t)
    constituent=average_Ehrhart_quasi_polynomial[0]#constituent of \tilde{L}^{gcd(n+1,period)}_{phi}(t)
    characteristic_polynomial_up_to_constant=shift_operator_act_on_polynomial(generalized_Eulerian_polynomial,constituent,n+1)#R(S^{n+1})constituent

    real_root_polynomial=shift_and_rot_polynomial(characteristic_polynomial_up_to_constant,int(Fraction(n*h,2)))#f(\sqrt{-1}*t+n*h/2)
    lb_nsol=Fourier_Budan(real_root_polynomial,-1000,1000)
    print("For the characteristic polynomial of the Linial arrangement \mathcal{A}",f"^[1,{n}],",sep="")
    if lb_nsol==len(characteristic_polynomial_up_to_constant)-1:
        print("All the roots of the characteristic polynomial have the same real part nh/2","\n")
    else:
        print("unknown")
