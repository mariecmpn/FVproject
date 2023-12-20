# Compilateur utilise
FC = gfortran-11


# Options de compilation
FFLAGS = -O0

# dossier pour les fichiers binaires
ODIR = ../bin

# chemin pour les fichiers executables
EXEDIR = ..

# Liste des fichiers binaires
OBJS1 = $(ODIR)/main.o $(ODIR)/schemas.o $(ODIR)/initialisation.o $(ODIR)/numerics.o
OBJS2 = $(ODIR)/main_2d.o $(ODIR)/schemas_2d.o $(ODIR)/initialisation_2d.o $(ODIR)/numerics.o

# Pour le modele scalaire
trafic_scalaire: $(OBJS1)
	$(FC) $(FFLAGS) $(OBJS1) -o $(EXEDIR)/trafic_scalaire

# Pour le modele systeme
trafic_system: $(OBJS2)
	$(FC) $(FFLAGS) $(OBJS2) -o $(EXEDIR)/trafic_system


# regle clean pour supprimer les fichiers binaires et executable
clean:
	rm $(ODIR)/*.o $(EXEDIR)/trafic_scalaire $(EXEDIR)/trafic_system

# Regles construction des fichiers binaires
$(ODIR)/main.o: main.f90 schemas.f90 initialisation.f90 numerics.f90 
	$(FC) -c $(FFLAGS) main.f90 -o $(ODIR)/main.o
$(ODIR)/numerics.o: numerics.f90 
	$(FC) -c $(FFLAGS) numerics.f90 -o $(ODIR)/numerics.o
$(ODIR)/schemas.o: schemas.f90 initialisation.f90 numerics.f90 
	$(FC) -c $(FFLAGS) schemas.f90 -o $(ODIR)/schemas.o
$(ODIR)/initialisation.o: initialisation.f90 numerics.f90 
	$(FC) -c $(FFLAGS) initialisation.f90 -o $(ODIR)/initialisation.o
$(ODIR)/main_2d.o: main_2d.f90 schemas_2d.f90 initialisation_2d.f90 numerics.f90 
	$(FC) -c $(FFLAGS) main_2d.f90 -o $(ODIR)/main_2d.o
$(ODIR)/schemas_2d.o: schemas_2d.f90 initialisation_2d.f90 numerics.f90 
	$(FC) -c $(FFLAGS) schemas_2d.f90 -o $(ODIR)/schemas_2d.o
$(ODIR)/initialisation_2d.o: initialisation_2d.f90 numerics.f90 
	$(FC) -c $(FFLAGS) initialisation_2d.f90 -o $(ODIR)/initialisation_2d.o
