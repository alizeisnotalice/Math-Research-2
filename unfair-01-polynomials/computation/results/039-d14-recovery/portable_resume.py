from pathlib import Path
import resume_gpu_factor_search_prefix as worker
worker.ORIGINAL=Path(__file__).resolve().parent/"frozen-worker/gpu_factor_search_prefix.py"
if __name__=="__main__": worker.main()
