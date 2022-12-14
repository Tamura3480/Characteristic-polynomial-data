### Abstract
Python program to compute the characteristic polynomial of the Linial arrangement of the exceptional root system and check the Postnikov-Stanley Linial arrangement conjecture. 

### Method
We transform the characteristic polynomial to the appropriate polynomial. Specifically, the characteristic polynomial is shifted by -nh/2 and then rotated 90 degrees around the origin. We count the number of real roots of the transformed characteristic polynomial using the Fourier Budan Theorem. If all the roots of the transformed characteristic polynomial are real, then all the roots of the characteristic polynomial have the same real part nh/2. 

### Usage
We use Python 3.8.5. 

The following two Python programs are available.    
・"compute_the_characteristic_polynomial.py"     
・"check_the_Postnikov_Stanley_Linial_arrangement_conjecture.py"

When we use them, it is necessary to assign "E_6","E_7","E_8" or "F_4" to the variable "phi" in each program. the variable "phi" means root system. 

When we determine the variable "phi" and run "compute_the_characteristic_polynomial.py", the characteristic polynomial of the Linial arrangement of "phi" is output to the file "characteristic_polynomial/characteristic_polynomial_{phi}.txt".

When we determine the variable "phi" and run "check_the_Postnikov_Stanley_Linial_arrangement_conjecture.py", we can check the conjecture for the Linial arrangement of "phi". If the conjecture holds, then the sentence "All the roots of the characteristic polynomial have the same real part nh/2" appears.

### Python3 program
* compute_the_characteristic_polynomial.py
* check_the_Postnikov_Stanley_Linial_arrangement_conjecture.py

### Folder
* constants_of_the_Ehrhart_quasi_polynomials.      
W:=Weyl group.      
f:=the index of connection.      
|W|/f is recorded in the files in this folder.      

* Ehrhart_quasi_polynomial_data.    
  The coefficient of the Ehrhart quasi polynomial, the period, and the number of divisor of the period are recorded in the files in this folder.

* Eulerian_polynomial_data.    
  The coefficient of the classical Eulerian polynomial in the files in this folder.

* root_system_ci_data.     
  The values of c_0, c_1, ・・・, c_l are recorded in the files in this folder.
