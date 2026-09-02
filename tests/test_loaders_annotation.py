"""Tests for the gene annotation file loaders."""

import pytest

from pylocuszoom.loaders import (
    load_bed,
    load_ensembl_genes,
    load_gtf,
)


class TestBEDLoader:
    """Tests for BED file loader."""

    def test_load_bed_basic(self, bed_file):
        """Test basic BED file loading."""
        df = load_bed(bed_file)

        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert len(df) == 3

    def test_load_bed_chromosome_cleaned(self, bed_file):
        """Test that chromosome prefix is removed."""
        df = load_bed(bed_file)

        # "chr1" should become "1"
        assert df["chr"].iloc[0] == "1"

    def test_load_bed12_format(self, tmp_path):
        """Test BED12 format with extra columns doesn't break."""
        # BED12 has: chr, start, end, name, score, strand, thickStart, thickEnd,
        #            itemRgb, blockCount, blockSizes, blockStarts
        content = """chr1\t1000000\t1020000\tGENE1\t100\t+\t1000500\t1019500\t0\t3\t100,200,300\t0,5000,19700
chr1\t1050000\t1080000\tGENE2\t200\t-\t1050000\t1080000\t0\t2\t500,600\t0,29400
"""
        filepath = tmp_path / "genes.bed12"
        filepath.write_text(content)

        df = load_bed(filepath)

        assert len(df) == 2
        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert df["gene_name"].iloc[0] == "GENE1"

    def test_load_bed_7_columns(self, tmp_path):
        """Test BED with 7 columns (more than 6)."""
        content = """chr1\t1000000\t1020000\tGENE1\t100\t+\textra_col
chr1\t1050000\t1080000\tGENE2\t200\t-\tmore_data
"""
        filepath = tmp_path / "genes.bed7"
        filepath.write_text(content)

        df = load_bed(filepath)

        assert len(df) == 2
        assert "gene_name" in df.columns
        assert df["gene_name"].iloc[0] == "GENE1"


class TestGTFLoader:
    """Tests for GTF file loader."""

    @pytest.fixture
    def gtf_file(self, tmp_path):
        """Create a temporary GTF file with gene_name before gene_id."""
        # Note: loader extracts first matching attribute (gene_name or gene_id)
        # so we put gene_name first to test gene_name extraction
        content = """##description: test GTF file
chr1\tENSEMBL\tgene\t1000000\t1020000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; gene_biotype "protein_coding";
chr1\tENSEMBL\texon\t1000000\t1005000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; exon_number 1;
chr1\tENSEMBL\texon\t1015000\t1020000\t.\t+\t.\tgene_name "BRCA1"; gene_id "ENSG00001"; exon_number 2;
chr1\tENSEMBL\tgene\t1050000\t1080000\t.\t-\t.\tgene_name "TP53"; gene_id "ENSG00002"; gene_biotype "protein_coding";
"""
        filepath = tmp_path / "genes.gtf"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def gtf_file_gene_id_only(self, tmp_path):
        """Create a GTF file with only gene_id (no gene_name)."""
        content = """##description: test GTF file
chr1\tENSEMBL\tgene\t1000000\t1020000\t.\t+\t.\tgene_id "ENSG00001"; gene_biotype "protein_coding";
"""
        filepath = tmp_path / "genes_id_only.gtf"
        filepath.write_text(content)
        return filepath

    def test_load_gtf_genes(self, gtf_file):
        """Test loading genes from GTF file."""
        df = load_gtf(gtf_file, feature_type="gene")

        assert len(df) == 2
        assert "chr" in df.columns
        assert "start" in df.columns
        assert "end" in df.columns
        assert "gene_name" in df.columns
        assert "strand" in df.columns

    def test_load_gtf_exons(self, gtf_file):
        """Test loading exons from GTF file."""
        df = load_gtf(gtf_file, feature_type="exon")

        assert len(df) == 2
        assert all(df["gene_name"] == "BRCA1")

    def test_load_gtf_chromosome_cleaned(self, gtf_file):
        """Test that chromosome prefix is removed."""
        df = load_gtf(gtf_file, feature_type="gene")

        # "chr1" should become "1"
        assert df["chr"].iloc[0] == "1"

    def test_load_gtf_extracts_gene_name(self, gtf_file):
        """Test that gene_name is extracted from attributes."""
        df = load_gtf(gtf_file, feature_type="gene")

        assert "BRCA1" in df["gene_name"].values
        assert "TP53" in df["gene_name"].values

    def test_load_gtf_preserves_strand(self, gtf_file):
        """Test that strand information is preserved."""
        df = load_gtf(gtf_file, feature_type="gene")

        brca1 = df[df["gene_name"] == "BRCA1"]
        tp53 = df[df["gene_name"] == "TP53"]
        assert brca1["strand"].iloc[0] == "+"
        assert tp53["strand"].iloc[0] == "-"

    def test_load_gtf_fallback_to_gene_id(self, gtf_file_gene_id_only):
        """Test that loader falls back to gene_id when gene_name missing."""
        df = load_gtf(gtf_file_gene_id_only, feature_type="gene")

        assert len(df) == 1
        assert df["gene_name"].iloc[0] == "ENSG00001"


class TestEnsemblGenesLoader:
    """Tests for Ensembl BioMart gene export loader."""

    @pytest.fixture
    def biomart_web_file(self, tmp_path):
        """Create a BioMart export using the web interface's display labels."""
        content = """Chromosome/scaffold name\tGene start (bp)\tGene end (bp)\tGene name\tStrand
1\t1000000\t1020000\tBRCA1\t1
1\t1050000\t1080000\tTP53\t-1
"""
        filepath = tmp_path / "biomart_web.tsv"
        filepath.write_text(content)
        return filepath

    @pytest.fixture
    def biomart_attr_file(self, tmp_path):
        """Create a BioMart export using the biomaRt attribute labels."""
        content = """chromosome_name\tstart_position\tend_position\texternal_gene_name\tstrand
1\t1000000\t1020000\tBRCA1\t1
1\t1050000\t1080000\tTP53\t-1
"""
        filepath = tmp_path / "biomart_attributes.tsv"
        filepath.write_text(content)
        return filepath

    def test_load_ensembl_genes_web_export_labels(self, biomart_web_file):
        """Test that BioMart display labels map to standard column names."""
        df = load_ensembl_genes(biomart_web_file)

        assert len(df) == 2
        assert df["chr"].iloc[0] == 1
        assert df["start"].iloc[0] == 1000000
        assert df["end"].iloc[0] == 1020000
        assert df["gene_name"].iloc[0] == "BRCA1"

    def test_load_ensembl_genes_web_export_strand_symbols(self, biomart_web_file):
        """Test that display-label strand integers become +/- symbols."""
        df = load_ensembl_genes(biomart_web_file)

        assert df["strand"].iloc[0] == "+"  # Not the raw 1
        assert df["strand"].iloc[1] == "-"  # Not the raw -1

    def test_load_ensembl_genes_attribute_labels(self, biomart_attr_file):
        """Test that biomaRt attribute labels map to standard column names."""
        df = load_ensembl_genes(biomart_attr_file)

        assert len(df) == 2
        assert df["chr"].iloc[0] == 1
        assert df["start"].iloc[0] == 1000000
        assert df["end"].iloc[0] == 1020000
        assert df["gene_name"].iloc[0] == "BRCA1"

    def test_load_ensembl_genes_attribute_strand_symbols(self, biomart_attr_file):
        """Test that attribute-label strand integers become +/- symbols."""
        df = load_ensembl_genes(biomart_attr_file)

        assert df["strand"].iloc[0] == "+"  # Not the raw 1
        assert df["strand"].iloc[1] == "-"  # Not the raw -1
