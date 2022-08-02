from vpython import *

#Fuerzas centrales!
#Daiana Dominikow
#Final Fisica Teorica 1

G=6.67e-11      #constante universal
rt=6.378e6      #radio de la tierra
mt=5.972e24     #masa de la tierra
ml=7.348e22     #masa de la luna
mu= mt*ml/(mt+ml)   #masa reducida
rl=1.7371e6     #radio de la luna
d=384.4e6       #distancia de la luna a la tierra

#primero creo las esferas y les doy sus atributos
tierra = sphere(pos = vector(0,0,0), 
                radius=rt*10,    #multiplique 10 veces el radio para que se vea mejor
                color=color.cyan, 
                make_trail=True)
                
luna = sphere(pos=tierra.pos+vector(d,0,0),
              radius=rl*10,
              make_trail=True,
              trail_type="points",
              interval=10, retain=500)

#la velocidad radial de la luna
vl=sqrt(G*mt/d)

#momento inicial del sistema
tierra.p = vector(0,0,0) 
luna.p = ml*vector(0,vl,0)  #p=m*v

#defino mi tiempo
t = 0
dt = 360
mes=24*3600*30

#momento angular
l=mag(cross(tierra.pos,tierra.p)+cross(luna.pos,luna.p))

#seteo el gráfico de la energía en funcion a la posición
grafico = graph(xtitle='Distancia (m)', ytitle = 'Energia (J)')
fU = gcurve(color=color.green, dot=True) 

#vamos a simular 6 meses con el metodo de Cromer - Euler
while t<mes*6:
    rate(5000) #velocidad de la simulacion
    r = luna.pos - tierra.pos
    F = -G*mt*ml*norm(r)/mag(r)**2 
    luna.p = luna.p + F*dt
    tierra.p = tierra.p - F*dt
    luna.pos = luna.pos + luna.p*dt/ml
    tierra.pos = tierra.pos + tierra.p*dt/mt
    
    #para graficar la energía potencial efectiva del sistema
    Ug = -G*ml*mt/mag(r)     #potencial gravitatoria
    Uc = l**2 / (2*mu*mag(r)**2)     #potencial centrípeta
    Uef = Ug+Uc     #potencial efectiva
    fU.plot(mag(r),Uef)
    t = t + dt