'''
LocalColabFold interface
'''

from pathlib import Path
from typing import Union, List, Tuple
from functools import cached_property

from .base import AlphaFoldBase
import colabfold.batch as cf
cf.logger.setLevel("INFO")


class ColabFold(AlphaFoldBase):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._queries = None
    def mmseq2(self, output_dir=None) -> None:
        '''
        Run MMseqs2 on with server.

        :param output_dir: Output directory. If set this will override any previously set output directory.
        :return: None
        '''

        if output_dir is None:
            output_dir = self._output_dir

        query, _ = self._get_query_from_fasta(self._fasta_string)

        cf.run(queries=query, 
            result_dir=output_dir, 
            num_models=0,
            is_complex=self.is_complex,
            user_agent="colabfold/1.5.5"
            )
        print("finished cf running stuff") 
        self.set_msa(Path(output_dir) / f"{self._name}.a3m")

    def predict(self, output_dir: str = None, msa="8:16", num_seeds=128, num_recycles=1):
        '''
        Run Alphafold on Colabfold.

        :param output_dir: Output directory. If set this will override any previously set output directory.
        :param msa: MSA range, e.g. '8:16'
        :param num_seeds: Number of seeds. Each seed will yield 5 structures from 5 models.
        :param num_recycles: Number of recycles
        '''
        print("running prediction!")
        if self._msa is not None:
            print("msa stuff!")
            self._queries = self._get_query_from_msa(self._msa)
            print("[colabfold] Found MSA input.")
        else:
            print(f"fasta stuff")
            self._queries, _ = self._get_query_from_fasta(self._fasta_string)
            print(self._queries)
        try:
            max_seq, max_extra_seq = msa.split(":")
        except ValueError as e:
            raise ValueError("Invalid msa argument. Please provide a valid range, e.g. '8:16'") from e
        print(f"almost to prediction part")
        if output_dir is None:
            output_dir = self._output_dir
        print(max_seq)
        print(max_extra_seq)
        return cf.run(queries=self._queries, 
                      result_dir=output_dir, 
                      is_complex=self.is_complex,
                      num_seeds=num_seeds,
                      num_models=5,
                      num_recycles=num_recycles,
                      #user_agent="colabfold/1.5.5",
                      max_seq=int(max_seq),
                      max_extra_seq=int(max_extra_seq),
                      )
    
    def _get_query_from_msa(self, a3m_string: str):

        (seqs, _) = cf.parse_fasta(a3m_string)
        if len(seqs) == 0:
            raise ValueError(f"Input MSA file is empty")
        query_sequence = seqs[0]
        # Use a list so we can easily extend this to multiple msas later
        a3m_lines = [a3m_string]
        queries = [(self._name, query_sequence, a3m_lines)]

        return queries
    
    @cached_property
    def is_complex(self):
        _, is_complex = self._get_query_from_fasta(self._fasta_string)
        return is_complex

    def _get_query_from_fasta(self, fasta_string: str):
        '''
        Get a query list from a single sequence fasta
        '''

        (sequences, headers) = cf.parse_fasta(fasta_string)
        queries = []
        for sequence, header in zip(sequences, headers):
            sequence = sequence.upper()
            if sequence.count(":") == 0:
                # Single sequence
                queries.append((header, sequence, None))
                is_complex = False
            else:
                # Complex mode
                queries.append((header, sequence.split(":"), None))
                is_complex = True
        
        return queries, is_complex
