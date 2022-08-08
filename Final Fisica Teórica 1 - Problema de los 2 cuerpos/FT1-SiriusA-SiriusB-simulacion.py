from vpython import *

#Fuerzas centrales!
#Daiana Dominikow
#Final Fisica Teorica 1

G=6.67e-11      #constante universal
rA= 1.1903e6     #radio Sirius A
rB= 5844      #radio Sirius B
mA=3.978e30     #masa Sirius A
mB=2.025e30     #masa Sirius B
mu= mA*mB/(mA+mB)   #masa reducida
d= 3e9       #distancia entre las estrellas

#posición centros de masa
xA = (mB/(mA+mB))*d
xB = -(mA/(mA+mB))*d

#primero creo las esferas y les doy sus atributos
SiriusA = sphere(pos = vector(xA,0,0), 
                radius=rA*500,    #multiplique 10 veces el radio para que se vea mejor
                color=color.cyan, 
                make_trail=True)
                
SiriusB = sphere(pos=vector(xB,0,0),
              radius=rB*500,
              color=color.red,
              make_trail=True
              )

#Centro de masa!
# cm = (mA*SiriusA.pos+mB*SiriusB.pos)/(mA+mB)
# CM = sphere(pos=cm,
#                 radius=0.1e9)

#la velocidad de Sirius B
#vB = sqrt(G*mA**2/(d*(mA+mB)))
vB = sqrt(G*mA*mag(SiriusB.pos)/(d**2))

#momento inicial del sistema
SiriusA.p = vector(0,0,0) 
SiriusB.p = mB*vector(0,vB,0)  #p=m*v

#defino mi tiempo
t = 0
dt = 100
t_f=24*3600*30

#momento angular
l=mag(cross(SiriusA.pos,SiriusA.p)+cross(SiriusB.pos,SiriusB.p))

#seteo el gráfico de la energía en funcion a la posición
grafico = graph(xtitle='Distancia (m)', ytitle = 'Energia (J)')
fU = gcurve(color=color.green, dot=True) 

# Metodo de Cromer - Euler
while t<t_f:
    rate(50) #velocidad de la simulacion
    r = SiriusB.pos - SiriusA.pos
    F = -G*mA*mB*norm(r)/mag(r)**2 
    SiriusB.p = SiriusB.p + F*dt
    SiriusA.p = SiriusA.p - F*dt
    SiriusB.pos = SiriusB.pos + SiriusB.p*dt/mB
    SiriusA.pos = SiriusA.pos + SiriusA.p*dt/mA
    # cm =  (mA*SiriusA.pos+mB*SiriusB.pos)/(mA+mB)#actualizo centro de masa
    # CM.pos = cm
    #para graficar la energía potencial efectiva del sistema
    Ug = -G*mA*mB/mag(r)            #potencial gravitatoria
    Uc = l**2 / (2*mu*mag(r)**2)    #potencial centrípeta
    Uef = Ug+Uc                     #potencial efectiva
    fU.plot(mag(r),Uef)
    t = t + dt