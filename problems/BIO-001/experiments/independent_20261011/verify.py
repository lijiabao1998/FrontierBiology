"""Independent CSV metadata audit; never import author code or evaluate expression."""
import argparse,copy,csv,datetime as dt,gzip,hashlib,json,os,subprocess,sys
from collections import Counter
from pathlib import Path
from admission_guard import require_admission
from common import compare

PATHS=(('cells_good_coverage',),('single_perturbation_cells',),('dual_perturbations',),
 ('cells_per_single_perturbation','min'),('cells_per_single_perturbation','max'),('cells_per_single_perturbation','p10'),('cells_per_single_perturbation','p90'),
 ('smallest_5','PLK4'),('smallest_5','HES7'),('smallest_5','ZBTB10'),('smallest_5','NIT1'),('smallest_5','CITED1'),
 ('largest_5','KLF1'),('largest_5','BAK1'),('largest_5','CEBPE'),('largest_5','UBASH3B'),('largest_5','ETS2'))
CONTROL_NAMES=tuple('PLUS_ONE_'+'_'.join(p) for p in PATHS)
def write(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
def targets(identity):
 return tuple(sorted({t for t in identity.split('_') if t and not t.isdigit() and not t.startswith('NegCtrl')}))
def project(audit):
 return {'.'.join(p):leaf(audit,p) for p in PATHS}
def leaf(obj,path):
 for k in path:obj=obj[k]
 return obj
def legacy(here,out,case,reference,round_record):
 destination=out/'legacy'/case;destination.mkdir(parents=True,exist_ok=False)
 claimed=destination/'claimed.json';write(claimed,reference)
 command=[sys.executable,'-X','utf8',str(here/'reference/legacy/bio001_r1/independent_check.py'),'--result',str(claimed),'--output',str(destination/'checker-result.json')]
 remaining=1800-(dt.datetime.now(dt.timezone.utc)-dt.datetime.fromisoformat(round_record['budget']['preparation_started_at'])).total_seconds()
 if remaining<=0:raise RuntimeError('whole preparation/execution budget exhausted')
 start=dt.datetime.now(dt.timezone.utc)
 with (destination/'stdout.txt').open('wb') as so,(destination/'stderr.txt').open('wb') as se:
  r=subprocess.run(command,stdout=so,stderr=se,timeout=remaining,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONUTF8='1'))
 evidence={'command':command,'started_at':start.isoformat(),'finished_at':dt.datetime.now(dt.timezone.utc).isoformat(),'exit_code':r.returncode}
 write(destination/'command.json',evidence)
 answer=json.loads((destination/'checker-result.json').read_text(encoding='utf-8'))
 checks=answer.get('checks',{})
 if set(checks)!=set(('cells_total','class_counts','single_perturbations','cells_per_single_median','fold_sizes_deterministic','spot_identities')) or any(type(v) is not bool for v in checks.values()) or not ((r.returncode==0 and answer.get('verdict')=='INDEPENDENT_CHECK_MATCH' and all(checks.values())) or (r.returncode==1 and answer.get('verdict')=='INDEPENDENT_CHECK_MISMATCH' and not all(checks.values()))):raise RuntimeError('legacy infrastructure/exit/verdict failure; no inferred rejection: '+case)
 evidence['accepted']=r.returncode==0 and answer.get('verdict')=='INDEPENDENT_CHECK_MATCH' and len(checks)==6 and all(v is True for v in checks.values())
 evidence['verdict']=answer.get('verdict');evidence['result_path']=str((destination/'checker-result.json').relative_to(out))
 return evidence
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--round',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 here=Path(__file__).resolve().parent;record=require_admission(a.repo,a.round,'BIO-001',here/'PROTOCOL.json')
 protocol=json.loads((here/'PROTOCOL.json').read_text(encoding='utf-8'))
 if tuple(tuple(x['path']) for x in protocol['controls'])!=PATHS or tuple(x['name'] for x in protocol['controls'])!=CONTROL_NAMES or any(x['operation']!='add' or type(x['delta']) is not int or x['delta']!=1 for x in protocol['controls']):raise RuntimeError('fixed 17-control contract changed')
 data=here/'reference/legacy/lit_data/GSE133344_filtered_cell_identities.csv.gz'
 if hashlib.sha256(data.read_bytes()).hexdigest()!='daf30337e7f6f07096d57e0d81db784bef00c87bd1fc927f018792c2f7af81e4':raise RuntimeError('metadata bytes changed')
 spots={'SET_KLF1__SET_KLF1':('KLF1','SET'),'NegCtrl0_BAK1__NegCtrl0_BAK1':('BAK1',),'NegCtrl10_NegCtrl0':(),'TGFBR2_IGDCC3_2':('IGDCC3','TGFBR2')}
 if any(targets(k)!=v for k,v in spots.items()):raise RuntimeError('token semantics failed')
 cells=0;good=0;classes=Counter();singles=Counter();duals=set();identities=Counter();seen=set()
 with gzip.open(data,'rt',encoding='utf-8',newline='') as f:
  reader=csv.DictReader(f)
  if reader.fieldnames!=protocol['csv_header']:raise RuntimeError('CSV header changed')
  barcode=reader.fieldnames[0]
  for row in reader:
   if None in row or any(v is None for v in row.values()):raise RuntimeError('malformed CSV row')
   if not row[barcode] or row[barcode] in seen:raise RuntimeError('empty/duplicate cell barcode')
   seen.add(row[barcode]);cells+=1
   if row['good_coverage'] not in ('True','False'):raise RuntimeError('unknown coverage boolean')
   good+=row['good_coverage']=='True';identity=row['guide_identity'];identities[identity]+=1;genes=targets(identity)
   if not genes:classes['control']+=1
   elif len(genes)==1:classes['single']+=1;singles[genes[0]]+=1
   elif len(genes)==2:classes['dual']+=1;duals.add(genes)
   else:raise RuntimeError('unexpected higher-order identity')
 counts=sorted(singles.values());n=len(counts)
 values={'cells_good_coverage':good,'single_perturbation_cells':sum(singles.values()),'dual_perturbations':len(duals),
  'cells_per_single_perturbation':{'min':counts[0],'max':counts[-1],'p10':counts[n//10],'p90':counts[9*n//10]},
  'smallest_5':dict(sorted(singles.items(),key=lambda kv:(kv[1],kv[0]))[:5]),'largest_5':dict(sorted(singles.items(),key=lambda kv:(-kv[1],kv[0]))[:5])}
 claimed=json.loads((here/'reference/author-results.json').read_text(encoding='utf-8'));baseline=compare(project(values),project(claimed['norman_audit']),0)
 anchors={'cells_total':cells,'class_counts':dict(classes),'single_perturbations':n,'median':counts[n//2],'p10_index':n//10,'p90_index':9*n//10,'all_rows_retained':cells==sum(classes.values())}
 anchor_expected={'cells_total':111445,'class_counts':{'single':57831,'dual':41759,'control':11855},'single_perturbations':105,'median':495,'p10_index':10,'p90_index':94,'all_rows_retained':True}
 anchor_check=compare(anchor_expected,anchors,0)
 out=a.output.parent;write(out/'independent-metadata.json',{'anchors':anchors,'verified_17_fields':project(values),'single_counts':dict(sorted(singles.items())),'distinct_dual_target_sets':len(duals)})
 old_clean=legacy(here,out,'CLEAN',claimed,record);controls=[]
 for name,path in zip(CONTROL_NAMES,PATHS):
  mutant=copy.deepcopy(claimed);obj=mutant['norman_audit']
  for key in path[:-1]:obj=obj[key]
  if type(obj[path[-1]]) is not int:raise RuntimeError('noninteger baseline claim')
  obj[path[-1]]+=1;check=compare(project(values),project(mutant['norman_audit']),0)
  old=legacy(here,out,name,mutant,record)
  controls.append({'name':name,'path':['norman_audit',*path],'operation':'add1','independent_accepted':check['accepted'],'independent_rejected':check['accepted'] is False,'semantic_errors':check['errors'],'legacy':old})
 fixed=tuple(x['name'] for x in controls)==CONTROL_NAMES and len(controls)==17
 ok=baseline['accepted'] is True and anchor_check['accepted'] is True and old_clean['accepted'] is True and fixed and all(x['independent_rejected'] is True for x in controls) and all(x['legacy']['accepted'] is True for x in controls)
 legacy_gap=sum(x['legacy']['accepted'] is True for x in controls)
 result={'problem':'BIO-001','baseline':baseline,'anchors':anchors,'anchor_check':anchor_check,'controls':controls,'control_count':17,'legacy_clean':old_clean,'legacy_accepted_mutants':legacy_gap,'legacy_coverage_gap_confirmed':legacy_gap==17,'all_criteria_passed':ok,
 'verdict':'INDEPENDENT_METADATA_17_FIELDS_PASS_WITH_LEGACY_COVERAGE_GAP' if ok and legacy_gap==17 else 'INDEPENDENT_METADATA_17_FIELDS_PASS' if ok else 'FAIL_CLOSED',
 'limits':protocol['limits']}
 write(a.output,result);print(json.dumps({'verdict':result['verdict'],'fields':17,'controls_rejected':sum(x['independent_rejected'] for x in controls),'legacy_accepted_mutants':legacy_gap,'all_criteria_pass':ok}));return 0 if ok else 1
if __name__=='__main__':
 try:sys.exit(main())
 except Exception:
  import traceback;traceback.print_exc();sys.exit(2)
