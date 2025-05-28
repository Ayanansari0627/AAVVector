import argparse
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Constants
AAV_LIMIT_BP = 4700  # Maximum base pairs per AAV vector
DEFAULT_OVERLAP_BP = 50
SUPPORTED_METHODS = ['intein', 'trans-splicing', 'custom']

class Tokenizer:
    """
    Tokenizer for prime editing system components.
    Splits full sequence into pegRNA, RT, and Cas9 segments by lengths.
    """
    def __init__(self, peg_len: int, rt_len: int):
        self.peg_len = peg_len
        self.rt_len = rt_len

    def tokenize(self, full_seq: Seq) -> dict:
        if len(full_seq) < self.peg_len + self.rt_len:
            raise ValueError("Full sequence shorter than specified component lengths.")
        peg = full_seq[:self.peg_len]
        rt = full_seq[self.peg_len:self.peg_len + self.rt_len]
        cas9 = full_seq[self.peg_len + self.rt_len:]
        return {'pegRNA': peg, 'RT': rt, 'Cas9': cas9}

class Splitter:
    """
    Splits components into two AAV vectors with optional overlap/linker.
    """
    def __init__(self, method: str = 'intein', overlap: int = DEFAULT_OVERLAP_BP):
        if method not in SUPPORTED_METHODS:
            raise ValueError(f"Unsupported method '{method}'. Choose from {SUPPORTED_METHODS}.")
        self.method = method
        self.overlap = overlap

    def split(self, parts: dict) -> dict:
        peg = parts['pegRNA']
        rt = parts['RT']
        cas9 = parts['Cas9']

        # Basic assignment: pegRNA+RT to vector A, Cas9 to vector B
        va = peg + rt
        vb = cas9

        # Apply overlap or linker for reconstitution
        if self.method == 'intein':
            overlap_seq = cas9[:self.overlap]
            va = va + overlap_seq
            vb = overlap_seq + vb
        elif self.method == 'trans-splicing':
            # Use generic splicing donor/acceptor
            donor = Seq('GGTGAG')  # example donor motif
            acceptor = Seq('AGGTGG')  # example acceptor motif
            va = va + donor + cas9[:self.overlap]
            vb = cas9[:self.overlap] + acceptor + vb
        elif self.method == 'custom':
            # Custom overlap equals first overlap bp of cas9
            custom_ol = cas9[:self.overlap]
            va = va + custom_ol
            vb = custom_ol + vb

        # Validate size constraints
        if len(va) > AAV_LIMIT_BP:
            raise ValueError(f"Vector A size {len(va)} bp exceeds limit of {AAV_LIMIT_BP} bp.")
        if len(vb) > AAV_LIMIT_BP:
            raise ValueError(f"Vector B size {len(vb)} bp exceeds limit of {AAV_LIMIT_BP} bp.")

        return {'VectorA': va, 'VectorB': vb}

class Reconstitutor:
    """
    Simulates reconstitution success probabilities based on method and overlap.
    """
    BASE_SUCCESS = {
        'intein': 0.85,
        'trans-splicing': 0.7,
        'custom': 0.5
    }

    def __init__(self, method: str, overlap: int):
        if method not in self.BASE_SUCCESS:
            raise ValueError(f"Method '{method}' not recognized.")
        self.method = method
        self.overlap = overlap

    def simulate(self) -> float:
        base = self.BASE_SUCCESS[self.method]
        # More overlap increases success up to a plateau
        adj = min(1.0, self.overlap / 100.0)
        noise = random.uniform(-0.05, 0.05)
        return round(base * adj + noise, 3)

def read_sequence(input_path: str) -> Seq:
    """
    Read DNA sequence from FASTA or plain text.
    """
    try:
        record = SeqIO.read(input_path, 'fasta')
        return record.seq
    except Exception:
        with open(input_path, 'r') as f:
            raw = f.read().strip().replace('\n', '')
            return Seq(raw)

def main():
    parser = argparse.ArgumentParser(
        description='AAV Vector Packaging Optimizer'
    )
    parser.add_argument('sequence', help='Path to input FASTA/text file with full sequence')
    parser.add_argument('--peg', type=int, default=200, help='Length of pegRNA segment (bp)')
    parser.add_argument('--rt', type=int, default=800, help='Length of reverse transcriptase segment (bp)')
    parser.add_argument('--method', choices=SUPPORTED_METHODS, default='intein', help='Reconstitution method')
    parser.add_argument('--overlap', type=int, default=DEFAULT_OVERLAP_BP, help='Overlap length (bp)')
    args = parser.parse_args()

    seq = read_sequence(args.sequence)
    print(f"✅ Sequence loaded: {len(seq)} bp")

    tok = Tokenizer(peg_len=args.peg, rt_len=args.rt)
    parts = tok.tokenize(seq)
    print(f"🔬 Tokenized lengths → pegRNA: {len(parts['pegRNA'])} bp, RT: {len(parts['RT'])} bp, Cas9: {len(parts['Cas9'])} bp")

    splitter = Splitter(method=args.method, overlap=args.overlap)
    vectors = splitter.split(parts)
    print(f"📦 Vector A: {len(vectors['VectorA'])} bp | Vector B: {len(vectors['VectorB'])} bp")

    if args.method == 'trans-splicing':
        print(f"🔁 Overlap region ({args.overlap} bp) added for trans-splicing")
    elif args.method == 'intein':
        print(f"🔗 Intein overlap ({args.overlap} bp) embedded")
    elif args.method == 'custom':
        print(f"🧩 Custom overlap region of {args.overlap} bp used")

    if len(vectors['VectorA']) <= AAV_LIMIT_BP and len(vectors['VectorB']) <= AAV_LIMIT_BP:
        print(f"✅ All segments fit within 4.7 kB packaging limit")
    else:
        print(f"❌ One or more vectors exceed the 4.7 kB packaging limit")

    rec = Reconstitutor(method=args.method, overlap=args.overlap)
    success = rec.simulate()
    print(f"📈 Simulated reconstitution success: ~{round(success * 100, 1)}%")


if __name__ == '__main__':
    import sys
    sys.argv = [
        'aav_optimizer.py',       # Dummy script name
        'sequence.txt',           # Your sequence file
        '--peg', '250',
        '--rt', '900',
        '--method', 'trans-splicing',
        '--overlap', '75'
    ]
    main()

