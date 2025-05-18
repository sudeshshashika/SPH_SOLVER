**************README****************
Here sloshing tank action is simulated using SPH method based on Matthias Müller's publication. Further, three layers of static particles are used to define the wall 
boundary action. Therefore, the first layer close to the fluid domin is placed 0.05*3 = 0.15 distance away from the fluid domain. The distance 0.15 m is also the smoothing 
length used in the files. 

There are two main cases studied here. change the input paths and output paths accordingly in the sph.cu file.
[1]. TANK3D --> A fluid voulume of 0.7 * 0.8 * 1.0 is represented using particles where distance between particles are equally spaced with  0.05 m.
                tank len_x = 0.84 , len_y = 0.76, len_x = 1.84.
                To visualize the breaking dam case comment out the sloshing action part inside the UPDATE loop. Refere the comment in kernel.h file


[2]. TANK3DDouble --> A fluid voulume of 0.7 * 0.8 * 1.0 is represented using particles where distance between particles are equally spaced with  0.05/2 = 0.025 m.
                tank len_x = 0.84 , len_y = 0.76, len_x = 1.84.
                According to my observtions this one have better results compared to 1st case.

[3]. OUTLINE--> This input can be used to generate wall outline. This will only generate the wall particle movements for sloshing action. Please be aware that to
                simulate this one for the exact same runtime for which ever the case you use above. therefore change the runtime parameter in OUTLINE.par file
                here plese comment out the CONTINITY AND MOMENTUM kernels when you generate the outline and also comment out the lines mentioned inside the UPDATE 
                kerenel. Please refere the note in the the kernel.h file

Above both the cases are tuned for the parameters specified in the .par files. Further, domainVolume is also used as 0.2. Hence, tuning is required  if you wish to 
change the parameters.

*************Wall boundary condition***********
Eventhough there is a seperation between one smoothing length (0.15 m) between wall and static particles, when fluid particles hit the right boundary (+Z) there are +0.01 to +0.02 penetration 
can be observed in the wall.
Further, when fluid becomes to the rest position there is a visual gap between some fluid particles while some fluid particles attached to the wall.

***********MAKE FILE****************
cmake version 3.16.3 in linux evironment is used when generating the makefile. 
$cmake -S. -B.

***********Reflective input files*********
If you wish to use reflective boundary condition, you can use these files as inputs. However, then you need to comment out sloshing action part and 
uncomment the boundaryCondition function inside UPDATE kernel. and also uncomment the  D.setboundaryLimits(0,6,0,6,0,6) and  D.setDampingFactors(0.15,0.15); in the  sph.cu

*********visualization paraview version 5.7.0 *****************
Run the outine case to get the moving wall ouline. And use ouline option to visualize the outer line.
import results files and use the flag as the scaler. Then select Edit the coloring and select the Enable opacity mapping for  suefaces to extract the fluid domain 


******Youtube Link  for the animation ************

https://youtu.be/uIfe4Ox2IMM
