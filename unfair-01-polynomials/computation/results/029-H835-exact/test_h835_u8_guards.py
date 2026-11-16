import os
"""Production-path semantic negative tests for the H8 u8 witness convention."""
import copy, hashlib, json, struct
from pathlib import Path
from tempfile import TemporaryDirectory
from verify_frontend import verify
from reverse_ball_kernel import balls_from_boxes
from verify_h835_certificate import galois_check, validate_row, validate_stream_metadata

ROOT=Path(__file__).resolve().parent
BALLS=ROOT/'H8-a3b5-root-balls.json'; LEDGER=Path(os.environ['H835_INPUT_DIR'])/'H8-a3b5-canonical-masks.u32le'
CERTIFICATE=ROOT/'certificate.json'; WITNESS=Path(os.environ['H835_INPUT_DIR'])/'witnesses.u8'
GALOIS=ROOT/'H8-mother-certificate.json'
def rejected(action):
 try: action()
 except (AssertionError,ValueError,KeyError,TypeError): return
 raise AssertionError('invalid production-path input accepted')
def synchronized(meta,witness):
 out=copy.deepcopy(meta);out['witness_bytes']=len(witness);out['witnesses_sha256']=hashlib.sha256(witness).hexdigest();return out
def main():
 front,boxes,canonical=verify(json.loads(BALLS.read_text()),LEDGER);meta,witness=json.loads(CERTIFICATE.read_text()),bytearray(WITNESS.read_bytes())
 validate_stream_metadata(meta,witness,canonical,BALLS,LEDGER);scale,factors=balls_from_boxes(boxes,96);words=[x[0] for x in struct.iter_unpack('<I',canonical)]
 standard=words.index(front['standard_mask']);nonstandard=next(i for i,w in enumerate(words) if w!=front['standard_mask'])
 leading=bytearray(witness);leading[2*nonstandard+1]=0;validate_stream_metadata(synchronized(meta,leading),leading,canonical,BALLS,LEDGER)
 rejected(lambda:validate_row(words[nonstandard],leading[2*nonstandard],leading[2*nonstandard+1],front['standard_mask'],factors,scale))
 nonstd255=bytearray(witness);nonstd255[2*nonstandard:2*nonstandard+2]=b'\xff\xff';validate_stream_metadata(synchronized(meta,nonstd255),nonstd255,canonical,BALLS,LEDGER)
 rejected(lambda:validate_row(words[nonstandard],255,255,front['standard_mask'],factors,scale))
 wrongstd=bytearray(witness);wrongstd[2*standard:2*standard+2]=b'\x00\x01';validate_stream_metadata(synchronized(meta,wrongstd),wrongstd,canonical,BALLS,LEDGER)
 rejected(lambda:validate_row(front['standard_mask'],0,1,front['standard_mask'],factors,scale))
 low=next(i for i,w in enumerate(words) if w.bit_count()==4);high=bytearray(witness);high[2*low]=0;high[2*low+1]=9;validate_stream_metadata(synchronized(meta,high),high,canonical,BALLS,LEDGER)
 rejected(lambda:validate_row(words[low],high[2*low],9,front['standard_mask'],factors,scale))
 truncated=witness[:-1];rejected(lambda:validate_stream_metadata(synchronized(meta,truncated),truncated,canonical,BALLS,LEDGER))
 with TemporaryDirectory(dir=Path(os.environ['H835_INPUT_DIR'])) as temp:
  bad=json.loads(GALOIS.read_text());bad['certificates'][0]['coefficients_ascending'][0]=True;p=Path(temp)/'wrong-h-bool.json';p.write_text(json.dumps(bad));rejected(lambda:galois_check(p))
 result={'status':'H835_production_u8_semantic_guards_passed','rejections':['synchronized_offset0_positive_leading','synchronized_nonstandard_255_255','synchronized_wrong_standard','synchronized_offset_exceeds_degree','synchronized_truncated_length','galois_H_bool']}
 Path(__file__).with_name('H835-u8-guard-tests.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
