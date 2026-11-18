from pathlib import Path
import json

D=Path(__file__).resolve().parent
a=json.loads((D/'macbook/12x5-B-p0.json').read_text())
b=json.loads((D/'transpose-p0.json').read_text())
assert (b['width'],b['height'],b['unsafe_edges'],b['prime_index'])==(5,12,'R',0)
fields=('peak_frequency_residues','unprotected_max_frequency_residues',
        'bad_second_frequency_residues','states_by_layer','processed_by_layer')
assert all(a[k]==b[k] for k in fields)
print('verified transposed top region modulo channel 0',sum(b['states_by_layer']))
