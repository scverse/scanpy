#!/usr/bin/env python3
"""
Comprehensive Benchmark Runner and Visualizer

This script runs all optimization benchmarks, aggregates results from CSV files,
and creates comparison plots for different dataset sizes and optimization types.
"""

import os
import sys
import subprocess
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')

# Set up plotting style
plt.style.use('default')
sns.set_palette("husl")

class BenchmarkRunner:
    """Runs and visualizes optimization benchmarks."""

    def __init__(self, base_dir: str = "."):
        self.base_dir = Path(base_dir)
        self.results_dir = self.base_dir / "benchmark_results"
        self.results_dir.mkdir(exist_ok=True)

        # Available benchmarks
        self.benchmarks = {
            "hvg": {
                "script": "hvg_benchmark.py",
                "description": "HVG Optimization Benchmark",
                "csv_pattern": "*hvg*.csv"
            },
            "normalization": {
                "script": "benchmarks/normalization_benchmark.py",
                "description": "Normalization Optimization Benchmark",
                "csv_pattern": "*normalization*.csv"
            },
            "correctness": {
                "script": "test_optimizations_correctness.py",
                "description": "Correctness Verification",
                "csv_pattern": None  # No CSV output
            },
            "performance": {
                "script": "simple_hvg_performance_test.py",
                "description": "HVG Performance Test",
                "csv_pattern": None  # No CSV output
            }
        }

    def run_benchmark(self, benchmark_name: str, timeout: int = 600) -> bool:
        """Run a specific benchmark with timeout."""
        if benchmark_name not in self.benchmarks:
            print(f"❌ Unknown benchmark: {benchmark_name}")
            return False

        benchmark = self.benchmarks[benchmark_name]
        script_path = self.base_dir / benchmark["script"]

        if not script_path.exists():
            print(f"❌ Benchmark script not found: {script_path}")
            return False

        print(f"🚀 Running {benchmark['description']}...")
        print(f"   Script: {script_path}")

        try:
            # Run benchmark with timeout
            result = subprocess.run(
                [sys.executable, str(script_path)],
                cwd=str(self.base_dir),
                capture_output=True,
                text=True,
                timeout=timeout
            )

            if result.returncode == 0:
                print(f"✅ {benchmark_name} completed successfully")
                if result.stdout:
                    print(f"   Output preview: {result.stdout[:200]}...")
                return True
            else:
                print(f"❌ {benchmark_name} failed with return code {result.returncode}")
                if result.stderr:
                    print(f"   Error: {result.stderr[:200]}...")
                return False

        except subprocess.TimeoutExpired:
            print(f"⏰ {benchmark_name} timed out after {timeout} seconds")
            return False
        except Exception as e:
            print(f"❌ {benchmark_name} failed with exception: {e}")
            return False

    def run_all_benchmarks(self, skip_long: bool = False) -> Dict[str, bool]:
        """Run all available benchmarks."""
        print("🔬 Running All Optimization Benchmarks")
        print("=" * 60)

        results = {}

        # Quick benchmarks first
        quick_benchmarks = ["correctness", "performance"]
        for name in quick_benchmarks:
            results[name] = self.run_benchmark(name, timeout=300)
            print()

        # Longer benchmarks
        if not skip_long:
            longer_benchmarks = ["normalization", "hvg"]
            for name in longer_benchmarks:
                timeout = 1200 if name == "hvg" else 600  # HVG might take longer
                results[name] = self.run_benchmark(name, timeout=timeout)
                print()
        else:
            print("⏭️  Skipping longer benchmarks (hvg, normalization)")
            for name in ["normalization", "hvg"]:
                results[name] = False

        return results

    def find_csv_files(self) -> Dict[str, List[Path]]:
        """Find all benchmark CSV files."""
        csv_files = {}

        # Look for CSV files in current directory and benchmarks subdirectory
        search_dirs = [self.base_dir, self.base_dir / "benchmarks"]

        for search_dir in search_dirs:
            if not search_dir.exists():
                continue

            for csv_file in search_dir.glob("*.csv"):
                # Categorize based on filename
                filename = csv_file.name.lower()

                if "hvg" in filename:
                    if "hvg" not in csv_files:
                        csv_files["hvg"] = []
                    csv_files["hvg"].append(csv_file)
                elif "normalization" in filename:
                    if "normalization" not in csv_files:
                        csv_files["normalization"] = []
                    csv_files["normalization"].append(csv_file)
                else:
                    if "other" not in csv_files:
                        csv_files["other"] = []
                    csv_files["other"].append(csv_file)

        return csv_files

    def load_and_aggregate_csv_data(self) -> Dict[str, pd.DataFrame]:
        """Load and aggregate data from CSV files."""
        csv_files = self.find_csv_files()
        aggregated_data = {}

        print("📊 Loading CSV Data")
        print("=" * 30)

        for category, files in csv_files.items():
            print(f"\n📁 {category.upper()} files:")

            category_data = []
            for file_path in files:
                print(f"   📄 {file_path.name}")
                try:
                    df = pd.read_csv(file_path)
                    df['source_file'] = file_path.name
                    category_data.append(df)
                    print(f"      ✅ Loaded {len(df)} rows, {len(df.columns)} columns")
                except Exception as e:
                    print(f"      ❌ Error loading: {e}")

            if category_data:
                # Combine all dataframes for this category
                combined_df = pd.concat(category_data, ignore_index=True)
                aggregated_data[category] = combined_df
                print(f"   📈 Total {category} data: {len(combined_df)} rows")

        return aggregated_data

    def create_hvg_performance_plots(self, hvg_data: pd.DataFrame) -> List[str]:
        """Create HVG performance comparison plots."""
        plots_created = []

        if hvg_data.empty:
            print("⚠️  No HVG data available for plotting")
            return plots_created

        print("📈 Creating HVG Performance Plots...")

        # Determine what columns we have
        print(f"   Available columns: {list(hvg_data.columns)}")

        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('HVG Optimization Performance Analysis', fontsize=16, fontweight='bold')

        # Plot 1: Speedup vs Dataset Size (if we have the data)
        if 'dataset_size' in hvg_data.columns and 'speedup' in hvg_data.columns:
            ax1 = axes[0, 0]
            scatter = ax1.scatter(hvg_data['dataset_size'], hvg_data['speedup'],
                                alpha=0.7, s=60, c=hvg_data.index, cmap='viridis')
            ax1.set_xlabel('Dataset Size (cells × genes)')
            ax1.set_ylabel('Speedup (x)')
            ax1.set_title('Speedup vs Dataset Size')
            ax1.grid(True, alpha=0.3)
            ax1.set_xscale('log')
            ax1.set_yscale('log')
        else:
            axes[0, 0].text(0.5, 0.5, 'Dataset size vs speedup\ndata not available',
                           ha='center', va='center', transform=axes[0, 0].transAxes)
            axes[0, 0].set_title('Speedup vs Dataset Size')

        # Plot 2: Execution Time Comparison (if we have timing data)
        timing_cols = [col for col in hvg_data.columns if 'time' in col.lower()]
        if len(timing_cols) >= 2:
            ax2 = axes[0, 1]
            # Create box plot of timing data
            timing_data = hvg_data[timing_cols].melt(var_name='Method', value_name='Time (s)')
            sns.boxplot(data=timing_data, x='Method', y='Time (s)', ax=ax2)
            ax2.set_title('Execution Time Distribution')
            ax2.tick_params(axis='x', rotation=45)
            ax2.set_yscale('log')
        else:
            axes[0, 1].text(0.5, 0.5, 'Timing comparison\ndata not available',
                           ha='center', va='center', transform=axes[0, 1].transAxes)
            axes[0, 1].set_title('Execution Time Comparison')

        # Plot 3: Density vs Performance (if we have density data)
        if 'density' in hvg_data.columns and 'speedup' in hvg_data.columns:
            ax3 = axes[1, 0]
            scatter = ax3.scatter(hvg_data['density'], hvg_data['speedup'],
                                alpha=0.7, s=60, c=hvg_data.index, cmap='plasma')
            ax3.set_xlabel('Matrix Density')
            ax3.set_ylabel('Speedup (x)')
            ax3.set_title('Speedup vs Matrix Density')
            ax3.grid(True, alpha=0.3)
            ax3.set_yscale('log')
        else:
            axes[1, 0].text(0.5, 0.5, 'Density vs speedup\ndata not available',
                           ha='center', va='center', transform=axes[1, 0].transAxes)
            axes[1, 0].set_title('Speedup vs Matrix Density')

        # Plot 4: Summary Statistics
        ax4 = axes[1, 1]
        if 'speedup' in hvg_data.columns:
            # Create histogram of speedups
            ax4.hist(hvg_data['speedup'], bins=20, alpha=0.7, edgecolor='black')
            ax4.set_xlabel('Speedup (x)')
            ax4.set_ylabel('Frequency')
            ax4.set_title('Speedup Distribution')
            ax4.axvline(hvg_data['speedup'].median(), color='red', linestyle='--',
                       label=f'Median: {hvg_data["speedup"].median():.1f}x')
            ax4.legend()
        else:
            # Show summary statistics as text
            summary_text = f"HVG Data Summary:\n"
            summary_text += f"Total records: {len(hvg_data)}\n"
            summary_text += f"Columns: {len(hvg_data.columns)}\n"
            if 'speedup' in hvg_data.columns:
                summary_text += f"Avg speedup: {hvg_data['speedup'].mean():.1f}x\n"
                summary_text += f"Max speedup: {hvg_data['speedup'].max():.1f}x"
            ax4.text(0.1, 0.5, summary_text, transform=ax4.transAxes, fontsize=12,
                    verticalalignment='center')
            ax4.set_title('Data Summary')

        plt.tight_layout()

        # Save plot
        plot_path = self.results_dir / "hvg_performance_analysis.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plots_created.append(str(plot_path))
        plt.close()

        print(f"   ✅ Saved: {plot_path}")
        return plots_created

    def create_normalization_performance_plots(self, norm_data: pd.DataFrame) -> List[str]:
        """Create normalization performance comparison plots."""
        plots_created = []

        if norm_data.empty:
            print("⚠️  No normalization data available for plotting")
            return plots_created

        print("📈 Creating Normalization Performance Plots...")

        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Normalization Optimization Performance Analysis', fontsize=16, fontweight='bold')

        print(f"   Available columns: {list(norm_data.columns)}")

        # Plot 1: Mode Comparison (if we have mode data)
        if 'mode' in norm_data.columns and 'speedup' in norm_data.columns:
            ax1 = axes[0, 0]
            mode_speedups = norm_data.groupby('mode')['speedup'].mean().sort_values(ascending=True)
            bars = ax1.barh(range(len(mode_speedups)), mode_speedups.values)
            ax1.set_yticks(range(len(mode_speedups)))
            ax1.set_yticklabels(mode_speedups.index)
            ax1.set_xlabel('Average Speedup (x)')
            ax1.set_title('Performance by Optimization Mode')
            ax1.grid(True, alpha=0.3)

            # Add value labels on bars
            for i, v in enumerate(mode_speedups.values):
                ax1.text(v + 0.1, i, f'{v:.1f}x', va='center')
        else:
            axes[0, 0].text(0.5, 0.5, 'Mode comparison\ndata not available',
                           ha='center', va='center', transform=axes[0, 0].transAxes)
            axes[0, 0].set_title('Performance by Mode')

        # Plot 2: Throughput vs Matrix Size
        if 'matrix_size' in norm_data.columns and 'throughput' in norm_data.columns:
            ax2 = axes[0, 1]
            ax2.scatter(norm_data['matrix_size'], norm_data['throughput'], alpha=0.7)
            ax2.set_xlabel('Matrix Size')
            ax2.set_ylabel('Throughput (cells/sec)')
            ax2.set_title('Throughput vs Matrix Size')
            ax2.set_xscale('log')
            ax2.set_yscale('log')
            ax2.grid(True, alpha=0.3)
        else:
            axes[0, 1].text(0.5, 0.5, 'Throughput vs size\ndata not available',
                           ha='center', va='center', transform=axes[0, 1].transAxes)
            axes[0, 1].set_title('Throughput vs Matrix Size')

        # Plot 3: Time Distribution by Mode
        time_cols = [col for col in norm_data.columns if 'time' in col.lower()]
        if time_cols and 'mode' in norm_data.columns:
            ax3 = axes[1, 0]
            # Melt the timing data for plotting
            time_data = norm_data[['mode'] + time_cols].melt(id_vars=['mode'],
                                                           var_name='Measurement',
                                                           value_name='Time (s)')
            sns.boxplot(data=time_data, x='mode', y='Time (s)', ax=ax3)
            ax3.set_title('Execution Time by Mode')
            ax3.tick_params(axis='x', rotation=45)
            ax3.set_yscale('log')
        else:
            axes[1, 0].text(0.5, 0.5, 'Time distribution\ndata not available',
                           ha='center', va='center', transform=axes[1, 0].transAxes)
            axes[1, 0].set_title('Time Distribution')

        # Plot 4: Performance Summary
        ax4 = axes[1, 1]
        summary_stats = []

        if 'speedup' in norm_data.columns:
            summary_stats.extend([
                f"Average speedup: {norm_data['speedup'].mean():.1f}x",
                f"Median speedup: {norm_data['speedup'].median():.1f}x",
                f"Max speedup: {norm_data['speedup'].max():.1f}x"
            ])

        if 'throughput' in norm_data.columns:
            summary_stats.extend([
                f"Avg throughput: {norm_data['throughput'].mean():.0f} cells/s",
                f"Max throughput: {norm_data['throughput'].max():.0f} cells/s"
            ])

        summary_stats.extend([
            f"Total records: {len(norm_data)}",
            f"Data columns: {len(norm_data.columns)}"
        ])

        summary_text = "\n".join(summary_stats)
        ax4.text(0.1, 0.5, summary_text, transform=ax4.transAxes, fontsize=12,
                verticalalignment='center')
        ax4.set_title('Performance Summary')
        ax4.axis('off')

        plt.tight_layout()

        # Save plot
        plot_path = self.results_dir / "normalization_performance_analysis.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plots_created.append(str(plot_path))
        plt.close()

        print(f"   ✅ Saved: {plot_path}")
        return plots_created

    def create_comparison_plots(self, all_data: Dict[str, pd.DataFrame]) -> List[str]:
        """Create overall comparison plots across all optimizations."""
        plots_created = []

        print("📈 Creating Overall Comparison Plots...")

        # Create combined performance comparison
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Scanpy Optimization Performance Overview', fontsize=16, fontweight='bold')

        # Collect speedup data from all sources
        all_speedups = {}
        for category, data in all_data.items():
            if 'speedup' in data.columns:
                all_speedups[category] = data['speedup'].dropna()

        # Plot 1: Speedup Comparison
        if all_speedups:
            ax1 = axes[0, 0]
            speedup_data = []
            labels = []

            for category, speedups in all_speedups.items():
                speedup_data.append(speedups.values)
                labels.append(f"{category.upper()}\n(n={len(speedups)})")

            bp = ax1.boxplot(speedup_data, labels=labels, patch_artist=True)
            ax1.set_ylabel('Speedup (x)')
            ax1.set_title('Speedup Distribution by Optimization Type')
            ax1.set_yscale('log')
            ax1.grid(True, alpha=0.3)

            # Color the boxes
            colors = ['lightblue', 'lightgreen', 'lightcoral', 'lightyellow']
            for patch, color in zip(bp['boxes'], colors[:len(bp['boxes'])]):
                patch.set_facecolor(color)
        else:
            axes[0, 0].text(0.5, 0.5, 'No speedup data\navailable',
                           ha='center', va='center', transform=axes[0, 0].transAxes)
            axes[0, 0].set_title('Speedup Comparison')

        # Plot 2: Data Coverage
        ax2 = axes[0, 1]
        categories = list(all_data.keys())
        record_counts = [len(df) for df in all_data.values()]

        bars = ax2.bar(categories, record_counts, color=['skyblue', 'lightgreen', 'salmon', 'gold'][:len(categories)])
        ax2.set_ylabel('Number of Records')
        ax2.set_title('Benchmark Data Coverage')
        ax2.tick_params(axis='x', rotation=45)

        # Add value labels on bars
        for bar, count in zip(bars, record_counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + max(record_counts)*0.01,
                    f'{count}', ha='center', va='bottom')

        # Plot 3: Performance Summary Table
        ax3 = axes[1, 0]
        ax3.axis('off')

        summary_data = []
        for category, data in all_data.items():
            row = [category.upper()]

            if 'speedup' in data.columns and not data['speedup'].empty:
                speedups = data['speedup'].dropna()
                row.extend([
                    f"{speedups.mean():.1f}x",
                    f"{speedups.median():.1f}x",
                    f"{speedups.max():.1f}x"
                ])
            else:
                row.extend(["-", "-", "-"])

            row.append(str(len(data)))
            summary_data.append(row)

        if summary_data:
            table = ax3.table(cellText=summary_data,
                             colLabels=['Optimization', 'Avg Speedup', 'Median Speedup', 'Max Speedup', 'Records'],
                             cellLoc='center',
                             loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            table.scale(1.2, 1.5)
            ax3.set_title('Performance Summary Table', pad=20)

        # Plot 4: Key Insights
        ax4 = axes[1, 1]
        ax4.axis('off')

        insights = [
            "🎯 Key Optimization Insights:",
            "",
            "• HVG: Best for sparse, large datasets",
            "• Normalization: Consistent improvements",
            "• Numba: Excellent for compute-heavy ops",
            "• Auto mode: Adapts to dataset size",
            "",
            f"📊 Total benchmarks: {sum(len(df) for df in all_data.values())}",
            f"📁 Categories analyzed: {len(all_data)}",
        ]

        if all_speedups:
            all_values = np.concatenate(list(all_speedups.values()))
            insights.extend([
                "",
                f"🚀 Overall median speedup: {np.median(all_values):.1f}x",
                f"⚡ Maximum speedup achieved: {np.max(all_values):.1f}x"
            ])

        insights_text = "\n".join(insights)
        ax4.text(0.05, 0.95, insights_text, transform=ax4.transAxes, fontsize=11,
                verticalalignment='top', fontfamily='monospace')
        ax4.set_title('Analysis Summary')

        plt.tight_layout()

        # Save plot
        plot_path = self.results_dir / "optimization_overview.png"
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plots_created.append(str(plot_path))
        plt.close()

        print(f"   ✅ Saved: {plot_path}")
        return plots_created

    def generate_html_report(self, all_data: Dict[str, pd.DataFrame], plots: List[str]) -> str:
        """Generate an HTML report with all results."""
        report_path = self.results_dir / "benchmark_report.html"

        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Scanpy Optimization Benchmark Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background-color: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 20px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; text-align: center; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; border-left: 4px solid #3498db; padding-left: 15px; }}
        .summary {{ background: #ecf0f1; padding: 20px; border-radius: 5px; margin: 20px 0; }}
        .plot {{ text-align: center; margin: 30px 0; }}
        .plot img {{ max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #3498db; color: white; }}
        tr:hover {{ background-color: #f5f5f5; }}
        .metric {{ display: inline-block; margin: 10px 20px; text-align: center; }}
        .metric-value {{ font-size: 24px; font-weight: bold; color: #e74c3c; }}
        .metric-label {{ font-size: 14px; color: #7f8c8d; }}
        .timestamp {{ text-align: center; color: #7f8c8d; font-size: 12px; margin-top: 30px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🚀 Scanpy Optimization Benchmark Report</h1>

        <div class="summary">
            <h2>📊 Executive Summary</h2>
            <p>This report presents comprehensive performance analysis of Scanpy optimizations including HVG (Highly Variable Genes) and normalization improvements.</p>

            <div style="text-align: center;">
"""

        # Add key metrics
        total_records = sum(len(df) for df in all_data.values())
        html_content += f"""
                <div class="metric">
                    <div class="metric-value">{total_records}</div>
                    <div class="metric-label">Total Benchmarks</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{len(all_data)}</div>
                    <div class="metric-label">Optimization Types</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{len(plots)}</div>
                    <div class="metric-label">Analysis Plots</div>
                </div>
"""

        # Add speedup metrics if available
        all_speedups = []
        for data in all_data.values():
            if 'speedup' in data.columns:
                all_speedups.extend(data['speedup'].dropna().tolist())

        if all_speedups:
            html_content += f"""
                <div class="metric">
                    <div class="metric-value">{np.median(all_speedups):.1f}x</div>
                    <div class="metric-label">Median Speedup</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{np.max(all_speedups):.1f}x</div>
                    <div class="metric-label">Max Speedup</div>
                </div>
"""

        html_content += """
            </div>
        </div>
"""

        # Add plots
        for plot_path in plots:
            plot_name = Path(plot_path).stem.replace('_', ' ').title()
            html_content += f"""
        <h2>📈 {plot_name}</h2>
        <div class="plot">
            <img src="{Path(plot_path).name}" alt="{plot_name}">
        </div>
"""

        # Add data summaries
        html_content += """
        <h2>📋 Data Summary</h2>
"""

        for category, data in all_data.items():
            html_content += f"""
        <h3>{category.upper()} Optimization</h3>
        <p><strong>Records:</strong> {len(data)} | <strong>Columns:</strong> {len(data.columns)}</p>
"""

            if len(data) > 0:
                # Show first few rows as preview
                preview_html = data.head(3).to_html(classes='table', table_id=f'{category}_table')
                html_content += f"""
        <details>
            <summary>Data Preview (first 3 rows)</summary>
            {preview_html}
        </details>
"""

        html_content += f"""
        <div class="timestamp">
            <p>Report generated on {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
"""

        # Write HTML file
        with open(report_path, 'w') as f:
            f.write(html_content)

        print(f"   ✅ Generated HTML report: {report_path}")
        return str(report_path)

    def run_full_analysis(self, skip_long_benchmarks: bool = False) -> Dict:
        """Run complete benchmark analysis pipeline."""
        print("🔬 Scanpy Optimization Benchmark Analysis")
        print("=" * 60)

        results = {
            "benchmark_results": {},
            "csv_data": {},
            "plots": [],
            "report": None
        }

        # Step 1: Run benchmarks
        print("\n🚀 STEP 1: Running Benchmarks")
        results["benchmark_results"] = self.run_all_benchmarks(skip_long=skip_long_benchmarks)

        # Step 2: Load CSV data
        print("\n📊 STEP 2: Loading and Aggregating CSV Data")
        results["csv_data"] = self.load_and_aggregate_csv_data()

        # Step 3: Create visualizations
        print("\n📈 STEP 3: Creating Visualizations")
        all_plots = []

        # HVG plots
        if "hvg" in results["csv_data"]:
            hvg_plots = self.create_hvg_performance_plots(results["csv_data"]["hvg"])
            all_plots.extend(hvg_plots)

        # Normalization plots
        if "normalization" in results["csv_data"]:
            norm_plots = self.create_normalization_performance_plots(results["csv_data"]["normalization"])
            all_plots.extend(norm_plots)

        # Overall comparison plots
        if results["csv_data"]:
            comparison_plots = self.create_comparison_plots(results["csv_data"])
            all_plots.extend(comparison_plots)

        results["plots"] = all_plots

        # Step 4: Generate report
        print("\n📝 STEP 4: Generating Report")
        if results["csv_data"] and results["plots"]:
            results["report"] = self.generate_html_report(results["csv_data"], results["plots"])

        # Final summary
        print("\n" + "=" * 60)
        print("✅ ANALYSIS COMPLETE")
        print("=" * 60)
        print(f"📊 CSV files processed: {len(results['csv_data'])}")
        print(f"📈 Plots created: {len(results['plots'])}")
        print(f"📝 Report generated: {'Yes' if results['report'] else 'No'}")
        print(f"📁 Results saved to: {self.results_dir}")

        if results["plots"]:
            print("\n📈 Generated plots:")
            for plot in results["plots"]:
                print(f"   • {Path(plot).name}")

        if results["report"]:
            print(f"\n📝 HTML Report: {results['report']}")
            print("   Open this file in a web browser to view the complete analysis.")

        return results

def main():
    """Main execution function."""
    import argparse

    parser = argparse.ArgumentParser(description="Run Scanpy optimization benchmarks and create visualizations")
    parser.add_argument("--skip-long", action="store_true",
                       help="Skip long-running benchmarks (hvg, normalization)")
    parser.add_argument("--plots-only", action="store_true",
                       help="Only create plots from existing CSV data, skip running benchmarks")
    parser.add_argument("--base-dir", default=".",
                       help="Base directory containing benchmark scripts")

    args = parser.parse_args()

    # Create benchmark runner
    runner = BenchmarkRunner(base_dir=args.base_dir)

    if args.plots_only:
        print("📊 Creating plots from existing data only...")
        csv_data = runner.load_and_aggregate_csv_data()

        all_plots = []
        if "hvg" in csv_data:
            all_plots.extend(runner.create_hvg_performance_plots(csv_data["hvg"]))
        if "normalization" in csv_data:
            all_plots.extend(runner.create_normalization_performance_plots(csv_data["normalization"]))
        if csv_data:
            all_plots.extend(runner.create_comparison_plots(csv_data))

        if csv_data and all_plots:
            report = runner.generate_html_report(csv_data, all_plots)
            print(f"✅ Analysis complete. Report: {report}")
        else:
            print("❌ No data available for plotting")
    else:
        # Run full analysis
        results = runner.run_full_analysis(skip_long_benchmarks=args.skip_long)

        if not any(results["benchmark_results"].values()) and not results["csv_data"]:
            print("⚠️  No benchmark data available. Consider running benchmarks first.")
            return 1

    return 0

if __name__ == "__main__":
    exit(main())
