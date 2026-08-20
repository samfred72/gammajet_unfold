TRIGGER=$1
SIM=${2:-pythia}
SYSTAG=${3:-nominal}
root -b -l -q "unfold.C(\"$TRIGGER\",\"$SIM\",\"$SYSTAG\")"
