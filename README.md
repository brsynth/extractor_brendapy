Extractor BrendaPy
==================

Extract the desired enzymatic data from [Brenda](https://www.brenda-enzymes.org/) 
using the [brendapy parser](https://github.com/matthiaskoenig/brendapy).
Here, we have chosen to store the datasets in json files.


Usage
-----

In the terminal

1. Create conda env
2. Activate conda env


list_parameter : _list_ _of_ _parameters_ _below_

path_files : _path_ _to_ _the_ _folder_ _containing_ _data_brenda.txt_ _and_
_file_ _result_ _json._

```
python -m extractor_brendapy --list_parameters ec unniprot ... --path_file_databrenda home/... --list_ec 1.1.1.1
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

```
git clone https://github.com/brsynth/extractor_brendapy.git
cd extractor_brendapy
conda env create -f environment.yaml -n <name>
```

Authors
-------
Nolwenn Paris
Joan Hérisson

Licence
-------