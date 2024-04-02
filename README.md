Extractor BrendaPy
==================

Extract the desired enzymatic data from [Brenda](https://www.brenda-enzymes.org/) 
using the [brendapy parser](https://github.com/matthiaskoenig/brendapy).
Here, we have chosen to store the datasets in json files.


Installation
------------

```
git clone https://github.com/brsynth/extractor_brendapy.git
cd extractor_brendapy
conda env create -f environment.yaml -n <name_env>
conda activate <name_env>
```

Usage
-----

The following 3 arguments are optional.

list_parameter : _list_ _of_ _parameters_ _below_

path_files : _path_ _to_ _the_ _folder_ _containing_ _data_brenda.txt_ _and_
_file_ _result_ _json._

list_ec : if this argument is not set, list of all known ec is generated

In the terminal

```
python -m extractor_brendapy --list_parameters [PARAM1 PARAM2 ...] --path_file_databrenda [PATH] --list_ec [PARAM1 PARAM2 ...]
```

List of possible parameters :

ec; uniprot; organism; ID; substrate; value; comment; units; refs; data; chebi; 
KM; KKM; KI; TN; IC50; ref; TS; SY; SU; ST; SP; SA; PU; NSP; MW; LO; GI; IN; 
CL; CF; AP; tissues; SN; RT; RN; RE


Authors
-------
* Nolwenn Paris
* Joan Hérisson

Licence
-------