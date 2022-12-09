# -*- coding: utf-8 -*-

#Python version 3.8.5
import math
from typing import List
from fractions import Fraction

#Compute the Eulerian polynomial A_l(t)
#input  :the root system phi (rank l)
#output :the Eulerian polynomial A_l(t) (rank l)
def generate_Eulerian_polynomial(phi :str) -> List[int]:
    if not(phi in {"E_6","E_7","E_8","F_4"}):
        raise ValueError("Only exception root system")
    with open(f"Eulerian_polynomial_data/Eulerian_polynomial_{phi}.txt") as polynomial_data:
        Eulerian_polynomial=list(map(int,polynomial_data.readline().split()))
    return Eulerian_polynomial


#Compute the Ehrhart quasi polynomial up to constant factor (|W|/f)*L_{phi}(t)  (f/|W| is the leading coefficient of L_{phi}(t))
#input  :the root system phi (rank l)
#output :the Ehrhart quasi polynomial of the fundamental alcove of phi up to constant factor (|W|/f)*L_{phi}(t)
def generate_Ehrhart_quasi_polynomial(phi :str) -> List[List[int]]:
    if not(phi in {"E_6","E_7","E_8","F_4"}):
        raise ValueError("Only exception root system")
    with open(f"Ehrhart_quasi_polynomial_data/Ehrhart_quasi_polynomial_{phi}.txt") as quasi_polynomial_data:
        period=int(quasi_polynomial_data.readline())
        ndivisor=int(quasi_polynomial_data.readline())
        Ehrhart_quasi_polynomial_data=[[0]]*period
        gcd=[0]*ndivisor
        for i in range(ndivisor):
            gcd[i]=int(quasi_polynomial_data.readline())
            Ehrhart_quasi_polynomial_data[gcd[i]-1]=list(map(int,quasi_polynomial_data.readline().split()))
    Ehrhart_quasi_polynomial=[[0]]*period
    for i in range(period):
        gcd=math.gcd(i+1,period)
        Ehrhart_quasi_polynomial[i]=Ehrhart_quasi_polynomial_data[gcd-1]
    return Ehrhart_quasi_polynomial

#Compute the quasi polynomial m*\tilde{f}^{k}(t) (m=period of f/gcd(k, period of f))
#input  :k and the the quasi_polynomial f(t)
#output :the quasi polynomial m*\tilde{f}^{k}(t)
def generate_average_quasi_polynomial(k :int, quasi_polynomial :List[List[int]]) -> List[List[int]]:
    period=len(quasi_polynomial)
    d=len(quasi_polynomial[0])
    if period%k!=0:
        raise ValueError
    average_quasi_polynomial=[]
    m=int(period/k)
    for i in range(k):
        sum_quasi_polynomial=[0]*d
        for j in range(m):
            sum_quasi_polynomial=list(map(lambda x,y:x+y, sum_quasi_polynomial,quasi_polynomial[(i+j*k)%period]))
        average_quasi_polynomial.append(sum_quasi_polynomial)
    return average_quasi_polynomial

#Compute the polynomial [c]_t=1+t+t^2+・・・+t^{c-1}
#input  :c
#output :[c]_t=1+t+t^2+・・・+t^{c-1}
def generate_cyc(c :int) -> List[int]:
    return list([1]*c)

#Compute the product of polynomials f(t)h(t)
#input  :the polynomial1 f(t) and the polynomial2 h(t)
#output :the polynomial f(t)h(t)
def product_polynomial(polynomial1 :List[int], polynomial2 :List[int]) -> List[int]:
    d1=len(polynomial1)
    d2=len(polynomial2)
    d=d1+d2-1
    polynomial1.extend([0]*(d-d1))
    polynomial2.extend([0]*(d-d2))
    new_polynomial=[0]*d
    for i in range(d):
        for j in range(i+1):
            new_polynomial[i]+=polynomial1[j]*polynomial2[i-j]
    return new_polynomial

#read data c_1,c_2,・・・,c_l of the root system phi
#input  :the root system phi (rank l)
#output :c_1,c_2,・・・,c_l of phi
def generate_c(phi :str) -> List[int]:
    with open(f"root_system_ci_data/root_system_ci_{phi}.txt") as c_data:
        c=list(map(int,c_data.readline().split()))
    return c

#Compute the generalized_Eulerian_polynomial R_{phi}(t)=[c_1]_t[c_2]_t・・・[c_l]_t A_l(t)
#input  :the root system phi
#output :the generalized Eulerian polynomial [c_1]_t[c_2]_t・・・[c_l]_t A_l(t)
def generate_generalized_Eulerian_polynomial(phi :str) -> List[int]:
    c=generate_c(phi)
    cyc=[1]
    for i in range(len(c)):
        cyc=product_polynomial(generate_cyc(c[i]),cyc)
    Eulerian_polynomial=generate_Eulerian_polynomial(phi)
    return product_polynomial(cyc,Eulerian_polynomial)

#Compute sum of polynomials f(t)+h(t)
#input  :the polynomial1 f(t) and the polynomial2 h(t)
#output :the polynomial f(t)+h(t)
def sum_polynomial(polynomial1 :List[int], polynomial2 :List[int]) -> List[int]:
    d1=len(polynomial1)
    d2=len(polynomial2)
    d=d1 if d2<=d1 else d2
    polynomial1.extend([0]*(d-d1))
    polynomial2.extend([0]*(d-d2))
    return list(map(lambda x,y:x+y,polynomial1,polynomial2))

#Compute the polynomial S^a t^index=(t-a)^index
#input  :a and index
#output :the polynomial (t-a)^index
def shift_power_act_on_t_power(a :int, index :int) -> List[int]:
    t_p=[0]*(index+1)
    a_power=1
    for i in range(0,index+1):
        comb_number=Fraction(math.factorial(index),math.factorial(i))
        comb_number=Fraction(comb_number,math.factorial(index-i))#long(factorial(index))//long(factorial(i))
        t_p[index-i]+=a_power*comb_number#long(long(a_power)*long(comb_number))
        a_power*=a#long(a)
    return t_p

#Compute the polynomial g(S^a)t^index
#input  :the shift operator g(S) ,index and a
#output :the polynomial g(S^a)t^index
def shift_operator_act_on_t_power(shift_operator :List[int] , index :int ,a :int=1) -> List[int]:
    d=len(shift_operator)
    new_polynomial=[0]*(index+1)
    for i in range(d):
        t_p=list(map(lambda t: shift_operator[i]*t,shift_power_act_on_t_power(-a*i,index)))
        new_polynomial=sum_polynomial(new_polynomial,t_p)
    return new_polynomial

#Compute the polynomial g(S^a)f(t)
#input  :the shift operator g(S) ,the polynomial f(t) and a
#output :the polynomial g(S^a)f(t)
def shift_operator_act_on_polynomial(shift_operator :List[int] ,polynomial :List[int] ,a :int=1) -> List[int]:
    d=len(polynomial)
    new_polynomial=[0]*d
    for i in range(d):
        t_p=list(map(lambda t:polynomial[i]*t,shift_operator_act_on_t_power(shift_operator,i,a)))
        new_polynomial=sum_polynomial(new_polynomial, t_p)
    return new_polynomial

#Compute the leading coefficient of m*(|W|/f)*\tilde{L_{phi}}^{n}(t) (p:=period of L_{phi} ,m:=p/gcd(n,p)
#input  :the root system phi , period and n
#output :m*|W|/f
def generate_constant(phi :str ,period :int ,n :int) -> int:
    with open(f"constant_of_Ehrhart_quasi_polynomial/constant_of_Ehrhart_{phi}.txt") as const_data:
        frac_W_f=int(const_data.readline().split()[0])
    m=Fraction(period,math.gcd(n,period))
    return m*frac_W_f

#Save the characteristic polynomial to characteristic_polynomial/characteristic_polynomial_{phi}.txt
#input  :the polynomial, the root system phi, n of of the Linial arrangement \mathcal{A}^[1,n], and mode of the file
#output :None
def save_polynomial_to_file(polynomial :List[int], phi :str, n :int, mode :str) -> None:
    polynomial.reverse()
    with open(f"characteristic_polynomial/characteristic_polynomial_{phi}", mode) as f:
        print("the characteristic polynomial of the Linial arrangement ",end="",file=f)
        print("\mathcal{A}",f"^[1,{n}] is as follow",sep="",file=f)
        i=len(polynomial)
        for coefficient in polynomial:
            if (coefficient>0) and len(polynomial)>i:
                print("+", end="", file=f)
            else:
                pass

            if coefficient==0:
                pass
            elif coefficient==1:
                print(f"t^{i-1}" ,end="", file=f)
            else:
                print(f"{coefficient}t^{i-1}" ,end="", file=f)
            i-=1
        print("",file=f)
        print("",file=f)

#main program
#Compute the characteristic polynomial of the Linial arrangement of the exceptional root system.
#When we determine the variable "phi" and run "compute_the_characteristic_polynomial.py", the characteristic polynomial of the Linial arrangement of "phi" is output to file "characteristic_polynomial/characteristic_polynomial_{phi}.txt".
#Only "E_6","E_7","E_8" or "F_4" can be selected as "phi"

phi="E_8"

rad_period={"E_6":6,"E_7":6,"E_8":30,"F_4":6}
arrangement_list=list(set([math.gcd(i,rad_period[phi])-1 for i in range(rad_period[phi])]))#list of n of Linial arrangement \mathcal{A}^[1,n]
with open(f"characteristic_polynomial/characteristic_polynomial_{phi}", "w") as f:
    pass

generalized_Eulerian_polynomial=generate_generalized_Eulerian_polynomial(phi)#R_{phi}(t)
Ehrhart_quasi_polynomial=generate_Ehrhart_quasi_polynomial(phi)#L_{phi}(t)
period=len(Ehrhart_quasi_polynomial)

for n in arrangement_list:
    average_Ehrhart_quasi_polynomial=generate_average_quasi_polynomial(math.gcd(n+1,period),Ehrhart_quasi_polynomial)#\tilde{L}^{gcd(n+1,period)}_{phi}(t)
    constituent=average_Ehrhart_quasi_polynomial[0]#constituent of \tilde{L}^{gcd(n+1,period)}_{phi}(t)
    characteristic_polynomial_up_to_constant=shift_operator_act_on_polynomial(generalized_Eulerian_polynomial,constituent,n+1)#R(S^{n+1})constituent

    constant=generate_constant(phi,period,n+1)#(period/gcd(period,n+1))*|W|/f
    characteristic_polynomial=list(map(lambda x:Fraction(x,constant), characteristic_polynomial_up_to_constant))

    save_polynomial_to_file(characteristic_polynomial,phi=phi ,n=n ,mode="a")
