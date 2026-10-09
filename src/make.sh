# Extra header directories (e.g. a RooUnfold checkout) from ROOT_INCLUDE_PATH
INC=$(echo "${ROOT_INCLUDE_PATH:-}" | tr ":" "\\n" | sed "/^$/d; s/^/-I/" | tr "\\n" " ")
rm -f unfolder.o drawer.o ana.o object.o pho_object.o jet_object.o treeuser.o unfold_utility.o insitu_utility.o purity_utility.o libgammajet_unfold.so
$(root-config --cxx) -c -fPIC -Wno-deprecated-declarations\
  drawer.cc \
  unfolder.cc \
  ana.cc \
  object.cc \
  pho_object.cc \
  jet_object.cc \
  treeuser.cc \
  unfold_utility.cc \
  insitu_utility.cc \
  purity_utility.cc \
  $INC `root-config --cflags`
echo ".o files made"
$(root-config --cxx) -shared -Wno-deprecated-declarations -o \
  libgammajet_unfold.so \
  drawer.o \
  unfolder.o \
  ana.o \
  object.o \
  pho_object.o \
  jet_object.o \
  treeuser.o \
  unfold_utility.o \
  insitu_utility.o \
  purity_utility.o \
  `root-config --libs`
echo "library made"
