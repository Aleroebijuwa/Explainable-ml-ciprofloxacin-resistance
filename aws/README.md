# Feature extraction on AWS EC2

`extract_features.py` downloads part of each isolate's sequencing reads from the European
Nucleotide Archive (ENA) and produces the three files used by `ecoli_ciprofloxacin_resistance.ipynb`:
`features.csv`, `qc.csv` and `mash_distances.csv`. It needs Linux bioinformatics tools, so it runs on
a rented Linux server rather than on Windows.

## 1. Launch a server
1. Sign in to the [AWS console](https://console.aws.amazon.com/) and choose the **Europe (London)
   eu-west-2** region (top right). London is close to ENA, so downloads are fast.
2. **Billing → Budgets → Create budget** (e.g. $10) so you get an email if costs go above it.
3. **EC2 → Launch instance**:
   - Name: `amr-features`
   - Image: **Ubuntu Server 24.04 LTS**
   - Instance type: **c6i.4xlarge** (16 CPUs, about $0.80/hour)
   - Key pair: **Create new key pair** → name `amr-key`, type RSA, format `.pem` → it downloads.
     Keep it; it is needed to copy the results back.
   - Storage: change 8 GiB to **40 GiB**
   - **Launch instance**

If launching fails with a *vCPU limit* error (common on new accounts), use **c6i.xlarge** (4 CPUs)
and `--workers 4`, or request a higher limit under **Service Quotas → EC2 → Running On-Demand
Standard instances**.

## 2. Run the extraction
Select the instance → **Connect** → **EC2 Instance Connect** → **Connect**. A terminal opens in the
browser. Paste:

```bash
git clone https://github.com/Aleroebijuwa/Explainable-ml-ciprofloxacin-resistance.git project
cd project
bash aws/setup.sh
nohup ~/venv/bin/python aws/extract_features.py --workers 16 > run.log 2>&1 &
tail -f run.log
```

`nohup ... &` keeps it running after the browser tab or laptop is closed. Press `Ctrl+C` to stop
watching the log (the run continues). To check progress later, connect again and run
`tail -f ~/project/run.log`. The script resumes where it stopped if it is started again.

It is finished when the log ends with `Written to /home/ubuntu/amr_cipro: features.csv ...`.

## 3. Copy the results to your laptop
In **PowerShell** on Windows, from the folder containing `amr-key.pem` (replace `<public-ip>` with
the instance's *Public IPv4 address* shown in EC2):

```powershell
scp -i amr-key.pem ubuntu@<public-ip>:~/amr_cipro/features.csv ubuntu@<public-ip>:~/amr_cipro/qc.csv ubuntu@<public-ip>:~/amr_cipro/mash_distances.csv .
```

If Windows complains that the key's permissions are too open, run once:
`icacls amr-key.pem /inheritance:r /grant:r "$($env:USERNAME):(R)"`

Put the three files in the project's `Data/` folder and run `ecoli_ciprofloxacin_resistance.ipynb`.

## 4. Terminate the server
**EC2 → Instances → select → Instance state → Terminate instance.** Stopping is not enough
(storage is still charged). Terminate once the three files are on your laptop.

## Options
| Option | Meaning |
|---|---|
| `--workers N` | isolates processed in parallel (default: number of CPUs) |
| `--n-isolates N` | process only a random subset, e.g. 20 for a quick test |
| `--no-retry-errors` | skip isolates that failed in an earlier run |
| `--combine-only` | only rebuild the three CSV files from finished isolates |
