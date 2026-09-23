"""Tests for Pydantic configuration classes.

Tests cover:
- RegionConfig: chrom >= 1, start >= 1, start < end, immutability
- ColumnConfig: sensible defaults, immutability
- DisplayConfig: sensible defaults, label_top_n >= 0, immutability
- LDConfig: lead_pos required when ld_reference_file provided, immutability
"""

import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, LDConfig, LDHeatmapInput
from pylocuszoom.exceptions import ValidationError
from pylocuszoom.plotter import LocusZoomPlotter
from tests.figure_probes import PROBES


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
        assert config.pos_col == "pos"
        assert config.p_col == "p_value"
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
        assert config.label_top_n is None
        assert config.auto_genes is None
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

    def test_with_defaults_fills_only_the_unset_fields(self):
        """None fields take the caller's defaults; set fields are kept."""
        from pylocuszoom.config import DisplayConfig

        unset = DisplayConfig().with_defaults(label_top_n=3, auto_genes=True)
        assert (unset.label_top_n, unset.auto_genes) == (3, True)
        chosen = DisplayConfig(label_top_n=0, auto_genes=False).with_defaults(
            label_top_n=3, auto_genes=True
        )
        assert (chosen.label_top_n, chosen.auto_genes) == (0, False)
        assert chosen.show_recombination is True

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
        with pytest.raises(ValidationError, match="needs a lead position"):
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

    def test_plot_config_defaults(self):
        """A region alone yields the documented column, display and LD defaults."""
        from pylocuszoom.config import PlotConfig, RegionConfig

        config = PlotConfig(region=RegionConfig(chrom=1, start=1000000, end=2000000))
        assert config.columns.pos_col == "pos"
        assert config.columns.p_col == "p_value"
        assert config.display.snp_labels is True
        assert config.ld.lead_pos is None
        assert config.panels.genes_df is None

    def test_plot_config_rejects_invalid_region(self):
        """The region rule fires when the composite is built."""
        from pylocuszoom.config import PlotConfig, RegionConfig

        with pytest.raises(ValidationError, match="start.*must be.*end"):
            PlotConfig(region=RegionConfig(chrom=1, start=2000, end=1000))


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
        """Each panel computes LD from its own fileset against its own lead."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000),
            n_panels=2,
            lead_positions=[1500, 1600],
            ld_reference_files=["/path/to/file1", "/path/to/file2"],
        )
        assert [(ld.lead_pos, ld.ld_reference_file) for ld in config.panel_lds()] == [
            (1500, "/path/to/file1"),
            (1600, "/path/to/file2"),
        ]

    def test_stacked_per_panel_ld_files_need_a_lead_per_panel(self):
        """The lead rule plot() applies holds for every stacked panel too."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        with pytest.raises(ValidationError, match="panel 1: .*needs a lead position"):
            StackedPlotConfig(
                region=RegionConfig(chrom=1, start=1000, end=2000),
                n_panels=2,
                ld_reference_files=["/path/to/file1", "/path/to/file2"],
            )

    def test_stacked_panel_ld_rejects_a_fileset_beside_a_broadcast_ld_col(self):
        """A per-panel fileset and a broadcast ld_col are the contradiction LDConfig bans."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        with pytest.raises(ValidationError, match="ld_col"):
            StackedPlotConfig(
                region=RegionConfig(chrom=1, start=1000, end=2000),
                n_panels=1,
                ld=LDConfig(ld_col="R2"),
                lead_positions=[1500],
                ld_reference_files=["/path/to/file1"],
            )

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

    def test_stacked_config_defaults_list_to_none(self):
        """List parameters should default to None, not empty lists."""
        from pylocuszoom.config import RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000, end=2000), n_panels=1
        )
        assert config.lead_positions is None
        assert config.panel_labels is None
        assert config.ld_reference_files is None

    def test_stacked_config_broadcast_ld_reference_file(self):
        """A broadcast ld_reference_file needs lead_positions, not ld.lead_pos.

        Bug fix: pyLocusZoom-vtf
        LD is computed per panel against that panel's lead, so the per-panel
        list satisfies the lead requirement.
        """
        from pylocuszoom.config import LDConfig, RegionConfig, StackedPlotConfig

        config = StackedPlotConfig(
            region=RegionConfig(chrom=1, start=1000000, end=2000000),
            n_panels=2,
            ld=LDConfig(ld_reference_file="/shared/plink_file"),
            lead_positions=[1500000, 1600000],
        )
        assert config.ld.ld_reference_file == "/shared/plink_file"
        assert config.ld.lead_pos is None
        assert config.lead_positions == [1500000, 1600000]

    def test_stacked_config_broadcast_ld_without_lead_positions_fails(self):
        """Broadcast ld_reference_file without lead_positions should still fail.

        If no lead_positions provided, there's no way to know which SNP to use
        as the LD reference for each panel.
        """
        from pylocuszoom.config import LDConfig, RegionConfig, StackedPlotConfig

        with pytest.raises(ValidationError, match="needs a lead position"):
            StackedPlotConfig(
                region=RegionConfig(chrom=1, start=1000000, end=2000000),
                n_panels=2,
                ld=LDConfig(ld_reference_file="/shared/plink_file"),
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


class TestRegionalOptionSurface:
    """The regional options are declared once, on the config models."""

    @staticmethod
    def _label_counts(fig):
        """The number of SNP labels drawn on each panel that carries any."""
        counts = [
            sum(1 for text in ax.texts if text.get_text())
            for ax in PROBES["matplotlib"].panels(fig)
        ]
        return [count for count in counts if count]

    @pytest.fixture
    def quiet(self):
        return DisplayConfig(show_recombination=False)

    def test_threshold_omitted_inherits_the_plotter(self, small_regional_gwas_df):
        plotter = LocusZoomPlotter(species="canine", genomewide_threshold=1e-5)
        fig = plotter.plot(
            small_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
        )
        assert PROBES["matplotlib"].hline_levels(fig) == pytest.approx([5.0])

    def test_threshold_none_draws_no_line(self, canine_plotter, small_regional_gwas_df):
        fig = canine_plotter.plot_stacked(
            [small_regional_gwas_df, small_regional_gwas_df],
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            significance_threshold=None,
        )
        assert PROBES["matplotlib"].panel_count(fig) == 2
        assert PROBES["matplotlib"].hline_levels(fig) == []

    def test_threshold_overrides_for_one_call(
        self, canine_plotter, small_regional_gwas_df
    ):
        fig = canine_plotter.plot(
            small_regional_gwas_df,
            chrom=1,
            start=1000000,
            end=2000000,
            display=DisplayConfig(show_recombination=False),
            significance_threshold=1e-6,
        )
        assert PROBES["matplotlib"].hline_levels(fig) == pytest.approx([6.0])
        assert canine_plotter.genomewide_threshold == 5e-8

    def test_label_top_n_takes_the_method_default(
        self, canine_plotter, small_regional_gwas_df, quiet
    ):
        single = self._label_counts(
            canine_plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=quiet,
            )
        )
        stacked = self._label_counts(
            canine_plotter.plot_stacked(
                [small_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                display=quiet,
            )
        )
        assert (single, stacked) == ([5], [3])

    def test_label_top_n_set_wins_on_both_methods(
        self, canine_plotter, small_regional_gwas_df
    ):
        display = DisplayConfig(show_recombination=False, label_top_n=2)
        single = self._label_counts(
            canine_plotter.plot(
                small_regional_gwas_df,
                chrom=1,
                start=1000000,
                end=2000000,
                display=display,
            )
        )
        stacked = self._label_counts(
            canine_plotter.plot_stacked(
                [small_regional_gwas_df],
                chrom=1,
                start=1000000,
                end=2000000,
                display=display,
            )
        )
        assert (single, stacked) == ([2], [2])

    def test_config_models_are_exported(self):
        import pylocuszoom

        for name in ("ColumnConfig", "DisplayConfig", "LDConfig", "PanelInputs"):
            assert name in pylocuszoom.__all__
            assert getattr(pylocuszoom, name) is getattr(
                __import__("pylocuszoom.config", fromlist=[name]), name
            )

    def test_from_kwargs_is_gone(self):
        from pylocuszoom.config import PlotConfig, StackedPlotConfig

        assert not hasattr(PlotConfig, "from_kwargs")
        assert not hasattr(StackedPlotConfig, "from_kwargs")


class TestPanelInputs:
    """The optional-panel inputs reject values the panels cannot draw."""

    def test_rejects_an_unknown_ld_heatmap_metric(self):
        with pytest.raises(ValidationError, match="metric"):
            LDHeatmapInput(matrix=pd.DataFrame([[1.0]]), snp_ids=["a"], metric="R2")

    @pytest.mark.parametrize(
        ("model", "field", "value"),
        [
            ("EqtlInput", "threshold", 5.0),
            ("EqtlInput", "threshold", 0.0),
            ("LDHeatmapInput", "height", -3),
            ("LDHeatmapInput", "height", 0),
        ],
    )
    def test_rejects_out_of_range_numbers(self, model, field, value):
        """A threshold outside (0, 1] or a non-positive height used to pass."""
        import pylocuszoom

        frame = pd.DataFrame([[1.0]])
        required = (
            {"data": frame}
            if model == "EqtlInput"
            else {"matrix": frame, "snp_ids": ["a"]}
        )
        with pytest.raises(ValidationError, match=field):
            getattr(pylocuszoom, model)(**required, **{field: value})

    def test_an_eqtl_option_needs_its_frame(self):
        """A gene filter with no eQTL frame is not a state PanelInputs can hold."""
        from pylocuszoom import EqtlInput

        with pytest.raises(ValidationError, match="data"):
            EqtlInput(gene="BRCA1")

    def test_frames_accept_a_spark_like_frame(self):
        """Every PanelInputs frame is collected through toPandas, as README promises."""
        from pylocuszoom import EqtlInput, PanelInputs

        class SparkLike:
            def __init__(self, frame):
                self.frame = frame

            def toPandas(self):
                return self.frame

        genes = pd.DataFrame({"chr": [1], "start": [1], "end": [2], "gene_name": ["G"]})
        eqtl = pd.DataFrame({"pos": [1], "p_value": [0.1]})

        panels = PanelInputs(
            genes_df=SparkLike(genes), eqtl=EqtlInput(data=SparkLike(eqtl))
        )

        assert panels.genes_df is genes
        assert panels.eqtl.data is eqtl

    def test_rejects_an_object_that_is_not_a_frame(self):
        from pylocuszoom import PanelInputs

        with pytest.raises(ValidationError, match="genes_df"):
            PanelInputs(genes_df=[1, 2, 3])


class TestUnknownFields:
    """A misspelt or removed option is an error, not a silently dropped keyword."""

    @pytest.mark.parametrize(
        ("model", "field"),
        [
            ("PanelInputs", "eqtl_df"),
            ("PanelInputs", "ld_heatmap_metric"),
            ("LDConfig", "ld_column"),
            ("DisplayConfig", "lable_top_n"),
            ("ColumnConfig", "chr_col"),
            ("GenomeWideStyle", "pallete"),
        ],
    )
    def test_unknown_field_raises_naming_it(self, model, field):
        import pylocuszoom

        with pytest.raises(ValidationError, match=field):
            getattr(pylocuszoom, model)(**{field: "x"})
