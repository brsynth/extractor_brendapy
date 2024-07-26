Extractor BrendaPy
==================

Extract the desired enzymatic data from [Brenda](https://www.brenda-enzymes.org/) 
using the [brendapy parser](https://github.com/matthiaskoenig/brendapy).
Here, we have chosen to store the datasets in json files.


TO DO
-----
- Complete writing the unit test functions for the RXN_CMP_file_creator.py file.
- Finish setting up the module for RXN_CMP from the args file of src2 and main from the terminal.


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

Please note:
When you want the parameters 'KM', 'KKM', 'KI', 'TN', 'IC50', you must add 
**substrate** to the list of requested elements.
Other elements that may be requested, by 5 parameters are : substrate, value,
comment, units, refs, data, chebi.

For parameters: "ref", "TS", "SY", "SU", "ST", "SP", "SA", "PU", "NSP", "MW", 
"LO", "GI", "IN", "CL", "CF", "AP" add **data**.
Other elements that may be requested, by the list of parameters above are :
data, refs, comment.


RXN AND CMP
===========

```
python -m RXM_CMP --path_file --input_file --file_RXM --file_RXM --file_CMP --mail --mdp 
```

--path_file :

--input_file : The initial JSON files, which include the reaction, substrates,
products, their stoichiometric coefficients, and whether the reaction is reversible.

--file_RXM : A dictionary where the values are lists of molecules.

--file_CMP : A dictionary where each key maps to a list of molecules, and each 
molecule in the list is associated with its SMILES representation.

--mail: mail for brenda
--mdp : password


Authors
-------
* Nolwenn Paris
* Joan Hérisson

Licence
-------