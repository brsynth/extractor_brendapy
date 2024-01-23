Brenda Enzyme
=============

Creating a json file, with a dataset from [Brenda](https://www.brenda-enzymes.org/)
using the [brendapy parser](https://github.com/matthiaskoenig/brendapy)

Usage
-----

In the terminal

1. Create conda env
2. Activate conda env

```
list_parameter=[ _list of parameters below_ ]
path_files= _path to the folder containing data_brenda.txt and file result json_
python -m brenda_enz_code list_parameter path_files
```

List of possible parameters :
* ec
* uniprot
* organism
* ID
* substrate
* value
* comment
* units
* refs
* data
* chebi
* KM
* KKM
* KI
* TN
* IC50
* ref
* TS
* SY
* SU
* ST
* SP
* SA
* PU
* NSP
* MW
* LO
* GI
* IN
* CL
* CF
* AP
* tissues
* SN
* RT
* RN
* RE

Installation
------------