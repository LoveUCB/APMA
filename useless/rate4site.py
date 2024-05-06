import subprocess

def run_rate4site(input_file,output_file):
    rate4site_command = f"rate4site -s {input_file} -o {output_file}"
    subprocess.run(rate4site_command, shell = True)


if __name__ == "__main__":
    run_rate4site("/home/wangjingran/APMA/data/query_msa.fasta","/home/wangjingran/APMA/data/score.txt")
