import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sp
import random
from matplotlib.colors import Normalize
from sympy import symbols, Eq, solve, cos ,sin , tan , asin , acos , pi , re
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import comb

#Donne_geo
diff_pot = 600
m = 9.1
x0 = 500
x1 = 900
x2 = 1500
c_c = (1.6*diff_pot)/((x1-x0)*m)
c_e = 0.05621
c_s = (1.6*diff_pot)/((x1-x0)*m)
pas = 2
dia = 10
R = dia/2
alpha = 0.13
limite = 599
E0 = 5
Resistance = 10**(8)
I_strip = diff_pot/Resistance
q = 1.6*10**(-7)
# Premiere etape

def plot_droite(m,b):
    # Générer des valeurs pour x sur une plage donnée
    x_values = range(-11,10 )  # Plage de valeurs de x à tracer

    # Calculer les valeurs correspondantes pour y en utilisant l'équation de la droite (y = mx + b)
    y_values = [(m * i) + (b) for i in x_values]

    # Tracer la droite
    plt.plot(x_values, y_values, label=f"y = {m}x + {b}")
    #plt.scatter(x, y, color='red', label="Point donné")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()

P_inf = 0.05
P_1e = 0.1
W = 60
Ee = 0
p = 1 # A revoir parce qu'il y a deux p

P_1r = 0.5
r = 1
Er = 30

e1 = 0.26
e2 = 2

r1 = 0.26
r2 = 2

t1 = 0.25
t2 = 0.8
t3 = 0.25
t4 = 1
gama_ts0 = 3.5
E_ts = 260

s = 1.54

# On definit les fonctions de Secondary electron yield 
def gama_e0(E):
    return P_inf + ( P_1e - P_inf )*math.exp( - ( ( abs( E - Ee ) / W )**p )/p )

def gama_r0(E):
    return P_1r*(1 - math.exp( - ( E/Er)**r  ))

def gama_e(E,teta):
    return gama_e0(E)*( 1 + e1*( 1  -  cos(teta)**e2  ))

def gama_r(E,teta):
    return gama_r0(E)*( 1 + r1*(  1- cos(teta)**r2  ))

def D(x,s):
    return (s*x) /( s - 1 + x**(s) )

def gama_teta(teta):
    return gama_ts0*( 1 + t1*(1-cos(teta)**t2))

def E_teta( teta ):
    return E_ts*( 1 + t3*(1-cos(teta)**t4))

def gama_ts(E, teta):
    return gama_teta(teta)*D((E/E_teta(teta)),s)

def gama_totale(E):
    return gama_ts(E,1.39) / (1 - gama_e(E,1.39) + gama_r(E,1.39))

# Calcul des Pn

# def P_prime(para,M,n):
#     return math.comb(M,n)*(para**(n))*((1-para)**(M-n))
def P_prime(para, M, n):
    return comb(M, n, exact=True) * (para**n) * ((1-para)**(M-n))

def P(gamae , gamar , para , M, n):
    if n == 0 :
        return ( 1 - gamae - gamar )*P_prime(para,M,n)
    if n == 1 :
        return ( 1 - gamae - gamar )*P_prime(para,M,n) + gamae + gamar
    else :
        return ( 1 - gamae - gamar )*P_prime(para,M,n) 
    

def Tirage(valeurs,probabilites):
    return  np.random.choice(valeurs, p=probabilites)

def Generate(E,teta):
    M = 20
    gamae = gama_e(E,teta)
    gamar = gama_r(E,teta)
    gamats = gama_ts(E,teta)
    gama_prime = gamats/ (1- gamae - gamar) 
    para = gama_prime/M
    probabilities = np.zeros(M+1)
    values = np.arange(M+1)
    for i in range(M+1):
        probabilities[i] = P(gamae , gamar , para , M, i)
    return Tirage(values , probabilities)


def Delta(a,b,c):
   return b*b-4*a*c

def Resolution (a,b,c):
    if Delta(a,b,c)<0 : 
        return False  # A revoir
    else :
        racine_delta = Delta(a,b,c)**0.5
        solution = (-b+racine_delta)/(2*a)
        return solution

def Resolution_bis (a,b,c):
    if Delta(a,b,c)<0 : 
        return False       #A revoir
    else :
        racine_delta = Delta(a,b,c)**0.5
        solution1 = (-b+racine_delta)/(2*a)
        solution2 = (-b-racine_delta)/(2*a)
        return [ solution1 , solution2 ]

def Check_iscts (iscts): # fonction qui distingue entre le cas simple ou le champs electrique est suppose constant par domaine ou variable
    if iscts == 1 :
        return True
    else :
        return False

def Rot( v , teta):
    matrix_3x3 = np.array( [ [ cos(teta) , -sin(teta) , 0], [ sin(teta) , cos(teta) , 0 ], [ 0 , 0 , 1 ] ] )
    return np.dot(matrix_3x3, v)

def generate_cosine_angle(num_samples):
    # Générer num_samples échantillons aléatoires entre 0 et 1 suivant une distribution uniforme
    u = np.random.uniform(size=num_samples)
    # Utiliser la transformation inverse pour obtenir les angles
    angles = np.arcsin(u)

    return angles

def Pho_ele(Energie_photon):  # Photo electron : A function that generates the initial electron matrix
    
    #angle = random.uniform(-pi/2 , pi/2 ) 
        
    norme_vitesse = (2*Energie_photon/m)**0.5 # velocity norm
    
    Mpe = np.zeros((3,3))
    #Mpe[1][0],Mpe[1][1] = norme_vitesse*cos(angle), norme_vitesse*sin(angle)
    Mpe[0][1] = 6
    Mpe[1][0],Mpe[1][1] = norme_vitesse, 0

    return Mpe

def premiere_arrivee ( Matrice_photo_electron  ): # A function taking the matrix from the previous fun and transporting the e to x0
    
    a = c_e/2    # soit au point du pore si le champs est suppose cst ou a x0 au cas contraire
    b = Matrice_photo_electron[1][0]
    c = -x0
    
    t = Resolution ( a , b , c ) # l'instant d'arrive du photo electron a la region I
    y = Matrice_photo_electron[1][1]*t + Matrice_photo_electron[0][1]
    vx = Matrice_photo_electron[1][0] + c_e*t
    
    Matrice_photo_electron[0][0], Matrice_photo_electron[0][1] = x0, y
    Matrice_photo_electron[1][0] = vx 
    Matrice_photo_electron[2][0], Matrice_photo_electron[2][1] = t , 0
    
    #print(Matrice_photo_electron)
    
    return Matrice_photo_electron

def Check_if_hit(Matrice_arrive): # Fonction qui verifie si l'electron rentre au pore
    
    # cette condition a pour but donner la partie entiere dans le cas negatif 
    
    if  (Matrice_arrive[0][1] - 1 ) < 0 :
        n = int ( (Matrice_arrive[0][1] - 1 ) / ( dia + pas ) ) - 1
    else :
        n = int ( (Matrice_arrive[0][1] - 1 ) / ( dia + pas ) )
    
    if  (Matrice_arrive[0][1] + 1 ) < 0 :
        m = int ( (Matrice_arrive[0][1] + 1 ) / ( dia + pas ) ) -1
    else :
        m = int ( (Matrice_arrive[0][1] + 1 ) / ( dia + pas ) )
    
    if n ==  (Matrice_arrive[0][1] - 1 ) / ( dia + pas ):
        return [False,n]
    elif n != m :
        return [False,n]
    elif n == m :
        return [True,n]

    
def Point_de_contact2(Mat,n,cts):
    Mat00, Mat01 , Mat02 , Mat10 , Mat11 , Mat12 =  Mat[0][0], Mat[0][1] , Mat[0][2] , Mat[1][0] , Mat[1][1] , Mat[1][2]
    x = symbols('x')
    polynome_4_deg_bis = ((cts**2)*(tan(alpha)**2)/4)*x**4 + (Mat10*cts*tan(alpha)**2 - Mat11*cts*tan(alpha))*x**3 + (Mat00*cts*tan(alpha)**2 - Mat01*cts*tan(alpha) + (Mat10**2)*tan(alpha)**2 - 2*Mat10*Mat11*tan(alpha) + Mat11**2 + Mat12**2 + cts*dia*n*tan(alpha) + cts*dia*tan(alpha)/2 + cts*n*pas*tan(alpha) + cts*pas*tan(alpha)/2 - cts*x0*tan(alpha)**2)*x**2 + (2*Mat00*Mat10*tan(alpha)**2 - 2*Mat00*Mat11*tan(alpha) - 2*Mat01*Mat10*tan(alpha) + 2*Mat01*Mat11 + 2*Mat02*Mat12 + 2*Mat10*dia*n*tan(alpha) + Mat10*dia*tan(alpha) + 2*Mat10*n*pas*tan(alpha) + Mat10*pas*tan(alpha) - 2*Mat10*x0*tan(alpha)**2 - 2*Mat11*dia*n - Mat11*dia - 2*Mat11*n*pas - Mat11*pas + 2*Mat11*x0*tan(alpha))*x + (Mat00**2)*tan(alpha)**2 - 2*Mat00*Mat01*tan(alpha) + 2*Mat00*dia*n*tan(alpha) + Mat00*dia*tan(alpha) + 2*Mat00*n*pas*tan(alpha) + Mat00*pas*tan(alpha) - 2*Mat00*x0*tan(alpha)**2 + Mat01**2 - 2*Mat01*dia*n - Mat01*dia - 2*Mat01*n*pas - Mat01*pas + 2*Mat01*x0*tan(alpha) + Mat02**2 - R**2 + (dia**2)*n**2 + (dia**2)*n + (dia**2)/4 + 2*dia*(n**2)*pas + 2*dia*n*pas - 2*dia*n*x0*tan(alpha) + dia*pas/2 - dia*x0*tan(alpha) + (n**2)*pas**2 + n*pas**2 - 2*n*pas*x0*tan(alpha) + (pas**2)/4 - pas*x0*tan(alpha) + (x0**2)*tan(alpha)**2
    racines = solve(polynome_4_deg_bis, x)
    solution = []
    for racine in racines :
        racine = re(racine)
        yf =  Mat01 + racine*Mat11
        xf = (cts/2)*racine**2 + Mat00 + racine*Mat10
        f = yf - tan(alpha)*xf + tan(alpha)*x0 - (pas + dia)*n - (pas + dia)/2
        zf = Mat02 + racine*Mat12
        final = (f**2 + zf**2)**0.5
        if abs(final - R) < 0.001 and racine > 0.001  :
            solution.append(racine) 
    #print('les solutions de contact dans ce passage est ', racines)
    if len(solution) == 0 :
        return True
    
    else :   
        t = min(solution)
        x = Mat[0][0] + Mat[1][0]*t +  cts*(t**2)/2 
        if x >= x1  :
            return False
        elif x <= x0 :
            return True            
        else  :
            return t

def generer_angles(n,teta1 , teta2 ): # Generating angles
    angles = np.zeros(n)
    for j in range(n):
        angle = random.uniform( teta1 , teta2 )
        angles[j] = angle
    return angles 


def emi_sec(Mat,n):
    
    # L'energie
    
    E = 0.5*m*( Mat[1][0]**2 + Mat[1][1]**2 + Mat[1][2]**2 )
    
    # Premiere etape determiner teta
    
    y_r = Mat[0][1] - tan(alpha)*Mat[0][0] + tan(alpha)*x0 - (pas + dia)*n - (pas + dia)/2
    z_r = Mat[0][2]
    #print('dans ce passage dans emi sec y et z sont', y , z)
    
    teta = math.atan2(z_r, y_r)

    if teta < 0:
        teta += 2 * math.pi
    
    Mat[0][1] , Mat[0][2] = R*cos(teta) + tan(alpha)*(Mat[0][0]-x0) +  (pas + dia)*n + (pas + dia)/2 , R*sin(teta)
    
    # deuxieme etape determiner er et eteta 
    
    e_r = np.array([0 , cos(teta) , sin(teta) ])
    e_o = np.array([0 , -sin(teta) , cos(teta) ])
    z = np.array([0 , 0 , 1 ])
    x = np.array([1 , 0 , 0 ])
    # Matrice de rotation
    
    normal = Rot(-e_r, alpha)
    x_rot  = Rot(x, alpha)
    mormal = np.cross(normal, x_rot)
    
    # produit scaliare 
    
    norme = ( Mat[1][0]**2 + Mat[1][1]**2 + Mat[1][2]**2 )**0.5
    vitesse = np.array([Mat[1][0]/norme , Mat[1][1]/norme ,Mat[1][2]/norme ])
    scalaire =  np.dot(normal , -vitesse ) 
    angle = acos(scalaire)
    if angle > 89*pi/180 :
        angle = 89*pi/180
    
    nombre_elec = Generate(E,angle)
    Resultat = []
    
    if nombre_elec != 0 :
        
        if nombre_elec == 1 :
            
            g_prime = gama_ts(E,angle)*(1 - gama_e(E,angle) - gama_r(E,angle) ) 
            
            P_e , P_r = gama_e(E,angle) , gama_r(E,angle)
            P_s = g_prime*math.exp(-g_prime)*( 1 - gama_e(E,angle) - gama_r(E,angle) )
            P_t = P_e + P_r + P_s
            
            P_e , P_r , P_s = P_e/P_t , P_r/P_t , P_s/P_t
            values = [0,1,2]
            probabilities = [P_e , P_r , P_s]
            process = Tirage(values , probabilities)
            
            if process == 0 :
                v_0 = ((2*E)/(m))**0.5
                c_n = cos(angle) 
                c_m =  np.dot(mormal , vitesse ) 
                c_x_rot =  np.dot(x_rot , vitesse )
                                
                u = c_n*normal + c_m*mormal + c_x_rot*x_rot
                
                # emission des elec secondaires
                
                M = np.zeros((3,3)) 
                
                M[0][0] , M[0][1] , M[0][2] = Mat[0][0] , Mat[0][1] , Mat[0][2]
                M[1][0] , M[1][1] , M[1][2] = v_0*u[0] , v_0*u[1] , v_0*u[2]
                M[2][0] , M[2][1] , M[2][2]  = Mat[2][0] , 0 , 0
                
                Resultat.append(M)
            
                return Resultat

            if process == 1 :
                Energie = random.uniform(E,E0)
                v_0 = ((2*Energie)/m)**0.5
            
            if process == 2 :
                v_0 = ((2*min(E,E0))/m)**0.5
            
            u0 = random.uniform(0.0174, 1)
            phi1 = np.arccos(u0)
            phi2 = random.uniform( 0 , 2*pi )
            
            u = cos(phi1)*normal + sin(phi2)*sin(phi1)*mormal + cos(phi2)*sin(phi1)*x_rot
            
            M = np.zeros((3,3)) 
            
            M[0][0] , M[0][1] , M[0][2] = Mat[0][0] , Mat[0][1] , Mat[0][2]
            M[1][0] , M[1][1] , M[1][2] = v_0*u[0] , v_0*u[1] , v_0*u[2]
            M[2][0] , M[2][1] , M[2][2]  = Mat[2][0] , 0 , 0
            
            Resultat.append(M)
            
            return Resultat
        
        else :
            
            v_0 = (((2*min((E/nombre_elec),E0)))/m)**0.5
            
            for i in range(nombre_elec):
                
                u0 = random.uniform(0.0174, 1)
                phi1 = np.arccos(u0)
                phi2 = random.uniform( 0 , 2*pi )
                
                u = cos(phi1)*normal + sin(phi2)*sin(phi1)*mormal + cos(phi2)*sin(phi1)*x_rot

                M = np.zeros((3,3)) 
                
                M[0][0] , M[0][1] , M[0][2] = Mat[0][0] , Mat[0][1] , Mat[0][2]
                M[1][0] , M[1][1] , M[1][2] = v_0*u[0] , v_0*u[1] , v_0*u[2]
                M[2][0] , M[2][1] , M[2][2]  = Mat[2][0] , 0 , 0
                
                Resultat.append(M)
            
            return Resultat
    else :
        return Resultat

def Recuperation( Mat ) : 
    
    # Fonction qui prend la matrice de l'electron a la sortie du pore
    # Donne la matrice de l'elctron au point de recuperation a l'anode
    
    a = c_s/2
    b = Mat[1][0]
    c = Mat[0][0] - x2 # x2 est x anode et x1 x fin mcp 
    
    t = Resolution( a , b , c )
    y2 = Mat[0][1] + Mat[1][1]*t
    vx = Mat[1][0] + c_s*t
    z2 = Mat[0][2] + Mat[1][2]*t
    M = np.zeros((3,3))

    M[0][0] , M[0][1] , M[0][2] = x2 ,  y2 , z2
    M[1][0] , M[1][1] , M[1][2] = vx , Mat[1][1] , Mat[1][2] 
    M[2][0] = Mat[2][0] + t
    
    return M

def Transporter1( Mat , t , cts) : 
    # Fonction qui prend la matrice de l'electron a la sortie du pore
    # Donne la matrice de l'elctron au point de recuperation a l'anode
    
    x = Mat[0][0] + Mat[1][0]*t + cts*(t**2)/2
    y2 = Mat[0][1] + Mat[1][1]*t
    z2 = Mat[0][2] + Mat[1][2]*t
    vx = Mat[1][0] + cts*t
    
    M = np.zeros((3,3))

    M[0][0] , M[0][1] , M[0][2] = x ,  y2 , z2
    M[1][0] , M[1][1] , M[1][2] = vx , Mat[1][1] , Mat[1][2]
    M[2][0] = Mat[2][0] + t
    M[2][1] = Mat[2][1] - t
    
    return M

def Transporter2( Mat , t , cts) : 
    # Fonction qui prend la matrice de l'electron a la sortie du pore
    # Donne la matrice de l'elctron au point de recuperation a l'anode
    
    x = Mat[0][0] + Mat[1][0]*t + cts*(t**2)/2
    y2 = Mat[0][1] + Mat[1][1]*t
    vx = Mat[1][0] + cts*t
    z2 = Mat[0][2] + Mat[1][2]*t
    
    M = np.zeros((3,3))

    M[0][0] , M[0][1] , M[0][2] = x ,  y2 , z2
    M[1][0] , M[1][1] , M[1][2] = vx , Mat[1][1] , Mat[1][2]
    M[2][0] = Mat[2][0] + t
    
    return M

def ajouter_element_trie(liste, element):
    if len(liste)>0:
        index = 0
        while index < len(liste) and liste[index][2][1] < element[2][1]:
            index += 1
        liste.insert(index, element)
    else :
        liste.append(element)


def enlever_premier_element(tableau):
    if len(tableau) > 0:
        return tableau[1:]
    else:
        return []

def enlever_elements_indices(tableau, indices):
    indices_tries = sorted(indices, reverse=True)
    for index in indices_tries:
        if 0 <= index < len(tableau):
            del tableau[index]
            #print(tableau)

def Erreur(Mat , n):
    y_r = (Mat[0][1] - tan(alpha)*Mat[0][0] + tan(alpha)*x0 - (pas + dia)*n - (pas + dia)/2)
    z_r = (Mat[0][2])
    if (y_r**2+z_r**2)**0.5 > R + 0.1 :
        return True
    else :
        return False
    
def Rearrangement( A1 , A2 , n , cts ):
    
    Emi = []
    N_Emi = []
    
    for M in A1 :
        time = Point_de_contact2(M, n , cts )
        if time ==   False :
            N_Emi.append(M)                    
        elif time ==   True :
            continue                   
        else :       
            M[2][1] =  time
            ajouter_element_trie(Emi , M )

    for M in A2 : 
        if M[0][0] > x1 : 
             N_Emi.append(M)
        else :
            time = Point_de_contact2(M, n , cts )
            if time ==   False :
                N_Emi.append(M)                    
            elif time ==   True :
                continue                   
            else :       
                M[2][1] =  time
                ajouter_element_trie(Emi , M )

    return Emi , N_Emi

def main( E , c ):
    cts = c
    Mat = Pho_ele ( E ) # la matrice du premier photo electron 
    Mat = premiere_arrivee(Mat) # la matrice au pore
    check = Check_if_hit ( Mat )
    Resultat_final = []
    
    if check[0] and -limite < Mat[0][1] < limite :                
        condition = 0
        Ensemble_emi = []
        for w in range(3):
            M = np.zeros((3,3))
            M[0][0] , M[0][1] , M[0][2] = Mat[0][0] - 0.1*w , Mat[0][1] , Mat[0][2]
            M[1][0] , M[1][1] , M[1][2] = Mat[1][0] , Mat[1][1] , Mat[1][2]
            M[2][0] = Mat[2][0] 
            time_1_hit = Point_de_contact2( M , check[1] , cts )
            M[2][1] = time_1_hit
            Ensemble_emi.append(M)
        Ensemble_non_emi = []
        variable = 0
        Etat = 0

        while condition == 0 :               
            if len(Ensemble_emi) != 0 :    
                # Block 1 pour le plot 
                I = 0
                lead_M = Ensemble_emi[0]
                temps = lead_M[2][1]
                instant = lead_M[2][0]
                lead_M = Transporter1( lead_M , temps , cts )
                Elec_secondaire = emi_sec(lead_M,check[1])
                Ensemble_emi = enlever_premier_element(Ensemble_emi) 
                
                if len(Ensemble_non_emi) != 0 : # On avance les non emissions avec l'increment de temps  
                    indice = []
                    for j in range(len(Ensemble_non_emi)):                               
                        M_passage = Transporter2( Ensemble_non_emi[j] , temps , cts )                  
                        if M_passage[0][0] >= x2 :      
                            M_passage = Recuperation( Ensemble_non_emi[j] )
                            Resultat_final.append(M_passage)
                            indice.append(j)
                        else :     
                            Ensemble_non_emi[j] = M_passage
                            if Ensemble_non_emi[j][0][0] < x1 :
                                I = I + q*Ensemble_non_emi[j][1][0]/(x1- x0) 
                    enlever_elements_indices(Ensemble_non_emi,indice)    

                for j in range(len(Ensemble_emi)): # On avance les prochaines emissions avec l'increment 
                    Ensemble_emi[j] = Transporter1( Ensemble_emi[j] , temps , cts )
                    I = I + q*Ensemble_emi[j][1][0]/(x1- x0)  

                for k in range(len(Elec_secondaire)):    # Ici on traite notre emission a l'instant donné
                    time = Point_de_contact2(Elec_secondaire[k], check[1] , cts )
                    if time ==   False :
                        Ensemble_non_emi.append(Elec_secondaire[k])                    
                    elif time ==   True :
                        continue                   
                    else :       
                        Elec_secondaire[k][2][1] =  time
                        ajouter_element_trie(Ensemble_emi , Elec_secondaire[k] )      

                if  Etat == 0 :

                    if I >= 0.05*I_strip :
                        instant_0 = instant
                        cts = cts*0.8
                        Ensemble_emi , Ensemble_non_emi = Rearrangement(Ensemble_emi  , Ensemble_non_emi , check[1] , cts )
                        Etat = 1
                else :

                    if instant >= instant_0 + 5 and  I < 0.05*I_strip :
                        instant_0 = instant
                        cts = cts/0.8
                        Ensemble_emi , Ensemble_non_emi = Rearrangement(Ensemble_emi  , Ensemble_non_emi , check[1] , cts )
                        if cts == c :
                            Etat = 0

                    if  instant >= instant_0 + 5 and  I > 0.05*I_strip :
                        cts = cts*0.8
                        instant_0 = instant
                        Ensemble_emi , Ensemble_non_emi = Rearrangement(Ensemble_emi  , Ensemble_non_emi , check[1] , cts )

            else :
                indice2 = []
                for j in range (len(Ensemble_non_emi) ) :
                    Ensemble_non_emi[j] = Recuperation( Ensemble_non_emi[j])
                    Resultat_final.append(Ensemble_non_emi[j])
                    indice2.append(j)
                enlever_elements_indices(Ensemble_non_emi,indice2)    

            if len(Ensemble_emi) == 0 and len(Ensemble_non_emi) == 0 :
                # Block 2
                condition= 1
    return Resultat_final

import pickle 

Rs = main(4.21,c_c)

nom_de_variable  = "Rs" 

with open("400µm600V_3.5_0.5_saturation_0.8_each_5ps.pkl", "wb") as fichier:     
    pickle.dump(Rs, fichier)


