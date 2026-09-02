"""Tests for Pydantic configuration classes.

Tests cover:
- RegionConfig: chrom >= 1, start >= 1, start < end, immutability
- ColumnConfig: sensible defaults, immutability
- DisplayConfig: sensible defaults, label_top_n >= 0, immutability
- LDConfig: lead_pos required when ld_reference_file provided, immutability
"""

import pytest
from pydantic import ValidationError


class TestRegionConfig:
    """Tests for RegionConfig validation and immutability."""

    def test_valid_region_creates_successfully(self):
        """Valid region parameters should create config."""
        from pylocuszoom.config import RegionConfig

        config = RegionConfig(chrom=1, start=1000, end=2000)
        assert config.chrom == 1
        assert config.start == 1000
        assert config.end == 2000

    def test_chrom_accepts_string_and_int(self):
        """Chromosome can be int or string (e.g., feline 'A1')."""
        from pylocuszoom.config import RegionConfig

        config_int = RegionConfig(chrom=1, start=1000, end=2000)
        assert config_int.chrom == 1

        config_str = RegionConfig(chrom="A1", start=1000, end=2000)
        assert config_str.chrom == "A1"

        config_x = RegionConfig(chrom="X", start=1000, end=2000)
        assert config_x.chrom == "X"

    def test_chrom_rejects_invalid_inputs(self):
        """Invalid chromosome values should raise ValidationError."""
        from pylocuszoom.config import RegionConfig

        with pytest.raises(ValidationError, match="must be >= 1"):
            RegionConfig(chrom=0, start=1000, end=2000)

        with pytest.raises(ValidationError, match="must be >= 1"):
            RegionConfig(chrom=-1, start=1000, end=2000)

        with pytest.raises(ValidationError, match="must not be empty"):
            RegionConfig(chrom="", start=1000, end=2000)

        with pytest.raises(ValidationError, match="must not be empty"):
            RegionConfig(chrom="  ", start=1000, end=2000)

    def test_start_must_be_less_than_end(self):
        """Start >= end should raise ValidationError."""
        from pylocuszoom.config import RegionConfig

        with pytest.raises(ValidationError, match="start.*must be.*end"):
            RegionConfig(chrom=1, start=2000, end=1000)

        with pytest.raises(ValidationError, match="start.*must be.*end"):
            RegionConfig(chrom=1, start=1000, end=1000)

    def test_start_is_one_based(self):
        """Coordinates are 1-based: start=0 is rejected, start=1 accepted."""
        from pylocuszoom.config import RegionConfig

        with pytest.raises(ValidationError, match="greater than or equal to 1"):
            RegionConfig(chrom=1, start=0, end=1000)

        config = RegionConfig(chrom=1, start=1, end=1000)
        assert config.start == 1

    def test_region_is_frozen(self):
        """Config should be immutable after creation."""
        from pylocuszoom.config import RegionConfig

        config = RegionConfig(chrom=1, start=1000, end=2000)
        with pytest.raises(ValidationError):
            config.start = 500


class TestColumnConfig:
    """Tests for ColumnConfig defaults and immutability."""

    def test_default_values_match_plotter_signature(self):
        """Default column names should match plotter.py defaults."""
        from pylocuszoom.config import ColumnConfig

        config = ColumnConfig()
        assert config.pos_col == "ps"
        assert config.p_col == "p_wald"
        assert config.rs_col == "rs"

    def test_custom_values_accepted(self):
        """Custom column names should be accepted."""
        from pylocuszoom.config import ColumnConfig

        config = ColumnConfig(pos_col="position", p_col="pvalue", rs_col="snp_id")
        assert config.pos_col == "position"
        assert config.p_col == "pvalue"
        assert config.rs_col == "snp_id"

    def test_column_config_is_frozen(self):
        """Config should be immutable after creation."""
        from pylocuszoom.config import ColumnConfig

        config = ColumnConfig()
        with pytest.raises(ValidationError):
            config.pos_col = "new_col"


class TestDisplayConfig:
    """Tests for DisplayConfig defaults, validation, and immutability."""

    def test_default_values_match_plotter_signature(self):
        """Default display settings should match plotter.py defaults."""
        from pylocuszoom.config import DisplayConfig

        config = DisplayConfig()
        assert config.snp_labels is True
        assert config.label_top_n == 5
        assert config.show_recombination is True
        assert config.figsize == (12.0, 8.0)

    def test_custom_values_accepted(self):
        """Custom display settings should be accepted."""
        from pylocuszoom.config import DisplayConfig

        config = DisplayConfig(
            snp_labels=False,
            label_top_n=10,
            show_recombination=False,
            figsize=(8.0, 6.0),
        )
        assert config.snp_labels is False
        assert config.label_top_n == 10
        assert config.show_recombination is False
        assert config.figsize == (8.0, 6.0)

    def test_label_top_n_must_be_non_negative(self):
        """label_top_n must be >= 0."""
        from pylocuszoom.config import DisplayConfig

        with pytest.raises(ValidationError, match="label_top_n"):
            DisplayConfig(label_top_n=-1)

    def test_label_top_n_zero_is_valid(self):
        """label_top_n of 0 is valid (means no labels)."""
        from pylocuszoom.config import DisplayConfig

        config = DisplayConfig(label_top_n=0)
        assert config.label_top_n == 0

    def test_display_config_is_frozen(self):
        """Config should be immutable after creation."""
        from pylocuszoom.config import DisplayConfig

        config = DisplayConfig()
        with pytest.raises(ValidationError):
            config.snp_labels = False


class TestLDConfig:
    """Tests for LDConfig validation and immutability."""

    def test_default_values(self):
        """Default LD config should have all None values."""
        from pylocuszoom.config import LDConfig

        config = LDConfig()
        assert config.lead_pos is None
        assert config.ld_reference_file is None
        assert config.ld_col is None

    def test_ld_reference_file_requires_lead_pos_in_plot_config(self):
        """ld_reference_file without lead_pos should raise in PlotConfig.

        Note: LDConfig itself doesn't validate because StackedPlotConfig needs
        to allow broadcast mode where lead_positions list is used instead.
        Validation happens at the composite config level.
        """
        from pylocuszoom.config import LDConfig, PlotConfig, RegionConfig

        # LDConfig alone doesn't raise (needed for broadcast mode)
        ld = LDConfig(ld_reference_file="/path/to/file")
        assert ld.ld_reference_file == "/path/to/file"
        assert ld.lead_pos is None

        # But PlotConfig should raise since single plots need lead_pos
        with pytest.raises(ValidationError, match="lead_pos.*required"):
            PlotConfig(
                region=RegionConfig(chrom=1, start=1000, end=2000),
                ld=ld,
            )

    def test_ld_reference_file_with_lead_pos_valid(self):
        """ld_reference_file with lead_pos should work."""
        from pylocuszoom.config import LDConfig

        config = LDConfig(lead_pos=1500, ld_reference_file="/path/to/file")
        assert config.lead_pos == 1500
        assert config.ld_reference_file == "/path/to/file"

    def test_ld_col_without_reference_file_valid(self):
        """Pre-computed LD column without reference file is valid."""
        from pylocuszoom.config import LDConfig

        config = LDConfig(ld_col="R2")
        assert config.ld_col == "R2"
        assert config.ld_reference_file is None

    def test_lead_pos_alone_valid(self):
        """lead_pos without reference file is valid (just highlight lead SNP)."""
        from pylocuszoom.config import LDConfig

        config = LDConfig(lead_pos=1500000)
        assert config.lead_pos == 1500000

    def test_ld_col_and_reference_file_mutually_exclusive(self):
        """Setting both ld_col and ld_reference_file should raise ValueError.

        ld_col means LD is pre-computed; ld_reference_file means compute LD.
        Both at once is contradictory.
        """
        from pylocuszoom.config import LDConfig

        with pytest.raises(
            ValidationError, match="Cannot specify both ld_col.*ld_reference_file"
        ):
            LDConfig(
                ld_col="R2",
                ld_reference_file="/path/to/file",
                lead_pos=1500,
            )

    def test_ld_col_alone_still_valid(self):
        """ld_col without ld_reference_file should work fine."""
        from pylocuszoom.config import LDConfig

        config = LDConfig(ld_col="R2")
        assert config.ld_col == "R2"
        assert config.ld_reference_file is None

    def test_ld_reference_file_alone_still_valid(self):
        """ld_reference_file without ld_col should work fine."""
        from pylocuszoom.config import LDConfig

        config = LDConfig(ld_reference_file="/path/to/file", lead_pos=1500)
        assert config.ld_reference_file == "/path/to/file"
        assert config.ld_col is None

    def test_ld_config_is_frozen(self):
        """Config should be immutable after creation."""
        from pylocuszoom.config import LDConfig

        config = LDConfig()
        with pytest.raises(ValidationError):
            config.lead_pos = 1000


class TestConfigIntegration:
    """Integration tests for config classes working together."""

    def test_all_configs_can_be_imported(self):
        """All config classes should be importable from config module."""
        from pylocuszoom.config import (
            ColumnConfig,
            DisplayConfig,
            LDConfig,
            RegionConfig,
        )

        assert RegionConfig is not None
        assert ColumnConfig is not None
        assert DisplayConfig is not None
        assert LDConfig is not None

    def test_configs_are_pydantic_models(self):
        """All configs should be Pydantic BaseModel subclasses."""
        from pydantic import BaseModel

        from pylocuszoom.config import (
            ColumnConfig,
            DisplayConfig,
            LDConfig,
            RegionConfig,
        )

        assert issubclass(RegionConfig, BaseModel)
        assert issubclass(ColumnConfig, BaseModel)
        assert issubclass(DisplayConfig, BaseModel)
        assert issubclass(LDConfig, BaseModel)

    def test_configs_support_model_dump(self):
        """Configs should support Pydantic v2 model_dump()."""
        from pylocuszoom.config import RegionConfig

        config = RegionConfig(chrom=1, start=1000, end=2000)
        dumped = config.model_dump()
        assert dumped == {"chrom": 1, "start": 1000, "end": 2000}

    def test_configs_support_model_copy(self):
        """Configs should support Pydantic v2 model_copy() for variations."""
        from pylocuszoom.config import DisplayConfig

        base = DisplayConfig()
        modified = base.model_copy(update={"figsize": (6.0, 4.0)})

        # Original unchanged
        assert base.figsize == (12.0, 8.0)
        # Copy has new value
        assert modified.figsize == (6.0, 4.0)


class TestPlotConfig:
    """Tests for PlotConfig composite class."""

    def test_plot_config_composes_all_configs(self):
        """PlotConfig should compose region, columns, display, and ld configs."""
        from pylocuszoom.config import (
            ColumnConfig,
            DisplayConfig,
            LDConfig,
            PlotConfig,
            RegionConfig,
        )

        config = PlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000),
        )
        # Defaults for other fields
        assert isinstance(config.region, RegionConfig)
        assert isinstance(config.columns, ColumnConfig)
        assert isinstance(config.display, DisplayConfig)
        assert isinstance(config.ld, LDConfig)

    def test_plot_config_with_all_nested_configs(self):
        """PlotConfig should accept all nested configs explicitly."""
        from pylocuszoom.config import (
            ColumnConfig,
            DisplayConfig,
            LDConfig,
            PlotConfig,
            RegionConfig,
        )

        config = PlotConfig(
            region=RegionConfig(chrom=5, start=5000, end=10000),
            columns=ColumnConfig(pos_col="position", p_col="pvalue"),
            display=DisplayConfig(snp_labels=False, label_top_n=10),
            ld=LDConfig(lead_pos=7500),
        )
        assert config.region.chrom == 5
        assert config.columns.pos_col == "position"
        assert config.display.snp_labels is False
        assert config.ld.lead_pos == 7500

    def test_plot_config_is_frozen(self):
        """PlotConfig should be immutable."""
        from pylocuszoom.config import PlotConfig, RegionConfig

        config = PlotConfig(region=RegionConfig(chrom=1, start=1000, end=2000))
        with pytest.raises(ValidationError):
            config.display = None

    def test_plot_config_from_kwargs_minimal(self):
        """from_kwargs should work with just region parameters."""
        from pylocuszoom.config import PlotConfig

        config = PlotConfig.from_kwargs(chrom=1, start=1000000, end=2000000)
        assert config.region.chrom == 1
        assert config.region.start == 1000000
        assert config.region.end == 2000000
        # Defaults should match plotter.py
        assert config.columns.pos_col == "ps"
        assert config.columns.p_col == "p_wald"
        assert config.display.snp_labels is True
        assert config.display.label_top_n == 5

    def test_plot_config_from_kwargs_with_ld_params(self):
        """from_kwargs should map LD parameters to LDConfig."""
        from pylocuszoom.config import PlotConfig

        config = PlotConfig.from_kwargs(
            chrom=1,
            start=1000000,
            end=2000000,
            lead_pos=1500000,
            ld_reference_file="/path/to/plink",
        )
        assert config.ld.lead_pos == 1500000
        assert config.ld.ld_reference_file == "/path/to/plink"

    def test_plot_config_from_kwargs_with_display_params(self):
        """from_kwargs should map display parameters to DisplayConfig."""
        from pylocuszoom.config import PlotConfig

        config = PlotConfig.from_kwargs(
            chrom=1,
            start=1000000,
            end=2000000,
            snp_labels=False,
            label_top_n=10,
            show_recombination=False,
            figsize=(8.0, 6.0),
        )
        assert config.display.snp_labels is False
        assert config.display.label_top_n == 10
        assert config.display.show_recombination is False
        assert config.display.figsize == (8.0, 6.0)

    def test_plot_config_from_kwargs_with_column_params(self):
        """from_kwargs should map column parameters to ColumnConfig."""
        from pylocuszoom.config import PlotConfig

        config = PlotConfig.from_kwargs(
            chrom=1,
            start=1000000,
            end=2000000,
            pos_col="position",
            p_col="pvalue",
            rs_col="snp_id",
        )
        assert config.columns.pos_col == "position"
        assert config.columns.p_col == "pvalue"
        assert config.columns.rs_col == "snp_id"

    def test_plot_config_from_kwargs_validates_on_construction(self):
        """from_kwargs should fail fast on invalid region."""
        from pylocuszoom.config import PlotConfig

        # Invalid region: start >= end
        with pytest.raises(ValidationError, match="start.*must be.*end"):
            PlotConfig.from_kwargs(chrom=1, start=2000, end=1000)

    def test_plot_config_from_kwargs_ld_col_param(self):
        """from_kwargs should accept ld_col for pre-computed LD."""
        from pylocuszoom.config import PlotConfig

        config = PlotConfig.from_kwargs(
            chrom=1, start=1000000, end=2000000, ld_col="R2"
        )
        assert config.ld.ld_col == "R2"


class TestStackedPlotConfig:
    """Tests for StackedPlotConfig with list-based parameters."""

    def test_stacked_config_has_list_parameters(self):
        """StackedPlotConfig should have lead_positions and panel_labels as lists."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000),
            n_panels=2,
            lead_positions=[1500, 1600],
            panel_labels=["Study A", "Study B"],
        )
        assert config.lead_positions == [1500, 1600]
        assert config.panel_labels == ["Study A", "Study B"]

    def test_stacked_config_ld_reference_files_list(self):
        """StackedPlotConfig should support multiple LD reference files."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000),
            n_panels=2,
            ld_reference_files=["/path/to/file1", "/path/to/file2"],
        )
        assert config.ld_reference_files == ["/path/to/file1", "/path/to/file2"]

    def test_stacked_config_single_ld_reference_file(self):
        """StackedPlotConfig should support single ld_reference_file for broadcast.

        Note: LDConfig requires lead_pos when ld_reference_file is provided.
        In practice, lead_positions list is used with stacked plots.
        """
        from pylocuszoom.config import LDConfig, RegionConfig, StackedPlotConfig

        # When using ld_reference_file in LDConfig, lead_pos is still required
        # This is because LD calculation needs a reference SNP
        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000),
            n_panels=2,
            ld=LDConfig(ld_reference_file="/shared/file", lead_pos=1500),
        )
        assert config.ld.ld_reference_file == "/shared/file"
        assert config.ld.lead_pos == 1500

    def test_stacked_config_is_frozen(self):
        """StackedPlotConfig should be immutable."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000), n_panels=1
        )
        with pytest.raises(ValidationError):
            config.lead_positions = [1500]

    def test_stacked_config_from_kwargs_minimal(self):
        """from_kwargs should work with just region parameters."""
        from pylocuszoom.config import StackedPlotConfig

        config = StackedPlotConfig.from_kwargs(
            n_panels=1, chrom=1, start=1000000, end=2000000
        )
        assert config.region.chrom == 1
        assert config.lead_positions is None
        assert config.panel_labels is None

    def test_stacked_config_from_kwargs_with_list_params(self):
        """from_kwargs should accept list parameters."""
        from pylocuszoom.config import StackedPlotConfig

        config = StackedPlotConfig.from_kwargs(
            n_panels=2,
            chrom=1,
            start=1000000,
            end=2000000,
            lead_positions=[1500000, 1600000],
            panel_labels=["Study A", "Study B"],
            ld_reference_files=["/path/a", "/path/b"],
        )
        assert config.lead_positions == [1500000, 1600000]
        assert config.panel_labels == ["Study A", "Study B"]
        assert config.ld_reference_files == ["/path/a", "/path/b"]

    def test_stacked_config_from_kwargs_validates_on_construction(self):
        """from_kwargs should fail fast on invalid region."""
        from pylocuszoom.config import StackedPlotConfig

        with pytest.raises(ValidationError, match="start.*must be.*end"):
            StackedPlotConfig.from_kwargs(n_panels=1, chrom=1, start=2000, end=1000)

    def test_stacked_config_from_kwargs_inherits_plot_config_params(self):
        """from_kwargs should accept all PlotConfig parameters too."""
        from pylocuszoom.config import StackedPlotConfig

        config = StackedPlotConfig.from_kwargs(
            n_panels=2,
            chrom=1,
            start=1000000,
            end=2000000,
            pos_col="position",
            p_col="pvalue",
            snp_labels=False,
            label_top_n=3,  # stacked default is 3
        )
        assert config.columns.pos_col == "position"
        assert config.columns.p_col == "pvalue"
        assert config.display.snp_labels is False
        assert config.display.label_top_n == 3

    def test_stacked_config_defaults_list_to_none(self):
        """List parameters should default to None, not empty lists."""
        from pylocuszoom.config import StackedPlotConfig

        config = StackedPlotConfig.from_kwargs(
            n_panels=1, chrom=1, start=1000, end=2000
        )
        assert config.lead_positions is None
        assert config.panel_labels is None
        assert config.ld_reference_files is None

    def test_stacked_config_from_kwargs_broadcast_ld_reference_file(self):
        """from_kwargs should accept broadcast ld_reference_file with lead_positions.

        Bug fix: pyLocusZoom-vtf
        When ld_reference_file is provided for broadcast and lead_positions list
        is provided, it should not require lead_pos in LDConfig.
        """
        from pylocuszoom.config import StackedPlotConfig

        # This should NOT raise - LD calculation will use lead_positions per panel
        config = StackedPlotConfig.from_kwargs(
            n_panels=2,
            chrom=1,
            start=1000000,
            end=2000000,
            ld_reference_file="/shared/plink_file",  # broadcast to all panels
            lead_positions=[1500000, 1600000],  # per-panel lead positions
        )
        assert config.ld.ld_reference_file == "/shared/plink_file"
        assert config.ld.lead_pos is None  # Not set at LDConfig level
        assert config.lead_positions == [1500000, 1600000]

    def test_stacked_config_from_kwargs_broadcast_ld_without_lead_positions_fails(self):
        """Broadcast ld_reference_file without lead_positions should still fail.

        If no lead_positions provided, there's no way to know which SNP to use
        as the LD reference for each panel.
        """
        from pydantic import ValidationError

        from pylocuszoom.config import StackedPlotConfig

        with pytest.raises(ValidationError, match="lead_pos.*required|lead_positions"):
            StackedPlotConfig.from_kwargs(
                n_panels=2,
                chrom=1,
                start=1000000,
                end=2000000,
                ld_reference_file="/shared/plink_file",  # broadcast
                # No lead_positions - should fail
            )


class TestColocConfig:
    """Tests for ColocConfig validation and immutability."""

    def test_default_values(self):
        """Test that default values are correct."""
        from pylocuszoom.config import ColocConfig

        config = ColocConfig()
        assert config.gwas_p_col == "p_gwas"
        assert config.eqtl_p_col == "p_eqtl"
        assert config.pos_col == "pos"
        assert config.rs_col == "rs"
        assert config.ld_col is None
        assert config.lead_snp is None
        assert config.gwas_threshold == 5e-8
        assert config.eqtl_threshold == 1e-5
        assert config.show_correlation is True
        assert config.color_by_effect is False
        assert config.gwas_effect_col is None
        assert config.eqtl_effect_col is None
        assert config.h4_posterior is None
        assert config.figsize == (8.0, 8.0)

    def test_threshold_validation(self):
        """Test that invalid thresholds raise ValidationError."""
        from pylocuszoom.config import ColocConfig

        # Threshold must be > 0
        with pytest.raises(ValidationError, match="gwas_threshold"):
            ColocConfig(gwas_threshold=0)

        with pytest.raises(ValidationError, match="gwas_threshold"):
            ColocConfig(gwas_threshold=-1e-8)

        # Threshold must be <= 1
        with pytest.raises(ValidationError, match="gwas_threshold"):
            ColocConfig(gwas_threshold=1.5)

        with pytest.raises(ValidationError, match="eqtl_threshold"):
            ColocConfig(eqtl_threshold=0)

    def test_h4_posterior_range(self):
        """Test that h4_posterior must be in [0, 1]."""
        from pylocuszoom.config import ColocConfig

        # Valid values at boundaries
        config_zero = ColocConfig(h4_posterior=0)
        assert config_zero.h4_posterior == 0

        config_one = ColocConfig(h4_posterior=1)
        assert config_one.h4_posterior == 1

        config_mid = ColocConfig(h4_posterior=0.95)
        assert config_mid.h4_posterior == 0.95

        # Invalid: < 0
        with pytest.raises(ValidationError, match="h4_posterior"):
            ColocConfig(h4_posterior=-0.1)

        # Invalid: > 1
        with pytest.raises(ValidationError, match="h4_posterior"):
            ColocConfig(h4_posterior=1.1)

    def test_effect_coloring_requires_columns(self):
        """Test color_by_effect=True without effect cols raises error."""
        from pylocuszoom.config import ColocConfig

        # Missing both columns
        with pytest.raises(
            ValidationError, match="color_by_effect.*requires.*gwas_effect_col"
        ):
            ColocConfig(color_by_effect=True)

        # Missing eqtl_effect_col
        with pytest.raises(
            ValidationError, match="color_by_effect.*requires.*eqtl_effect_col"
        ):
            ColocConfig(color_by_effect=True, gwas_effect_col="beta_gwas")

        # Missing gwas_effect_col
        with pytest.raises(
            ValidationError, match="color_by_effect.*requires.*gwas_effect_col"
        ):
            ColocConfig(color_by_effect=True, eqtl_effect_col="beta_eqtl")

        # Valid: both columns provided
        config = ColocConfig(
            color_by_effect=True,
            gwas_effect_col="beta_gwas",
            eqtl_effect_col="beta_eqtl",
        )
        assert config.color_by_effect is True
        assert config.gwas_effect_col == "beta_gwas"
        assert config.eqtl_effect_col == "beta_eqtl"

    def test_frozen_config(self):
        """Test that ColocConfig is immutable after creation."""
        from pylocuszoom.config import ColocConfig

        config = ColocConfig()
        with pytest.raises(ValidationError):
            config.gwas_p_col = "new_col"

    def test_custom_values_accepted(self):
        """Test that custom values are accepted."""
        from pylocuszoom.config import ColocConfig

        config = ColocConfig(
            gwas_p_col="pvalue",
            eqtl_p_col="pval_eqtl",
            pos_col="position",
            rs_col="snp_id",
            ld_col="r2",
            lead_snp="rs12345",
            gwas_threshold=1e-5,
            eqtl_threshold=1e-3,
            show_correlation=False,
            figsize=(10.0, 10.0),
        )
        assert config.gwas_p_col == "pvalue"
        assert config.eqtl_p_col == "pval_eqtl"
        assert config.pos_col == "position"
        assert config.rs_col == "snp_id"
        assert config.ld_col == "r2"
        assert config.lead_snp == "rs12345"
        assert config.gwas_threshold == 1e-5
        assert config.eqtl_threshold == 1e-3
        assert config.show_correlation is False
        assert config.figsize == (10.0, 10.0)
