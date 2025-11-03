"""
SPICE Netlist Graph Coloring Analyzer
Determines the minimum number of colors needed for conflict-free parallel execution of MOSFETs

Usage: python3 color.py netlist.cir
"""

import sys
import os
import re
from collections import defaultdict
from typing import Set, List, Tuple, Dict

class SPICEColoringAnalyzer:
    def __init__(self):
        self.mosfets = []
        self.conflict_graph = defaultdict(set)
        self.colors = {}
        self.num_colors = 0
    
    def parse_netlist(self, netlist_content: str) -> List[Dict]:
        """Parse SPICE netlist and extract MOSFET devices."""
        mosfets = []
        lines = netlist_content.strip().split('\n')
        
        # Join continued lines (lines starting with '+')
        processed_lines = []
        current_line = ""
        for line in lines:
            line = line.strip()
            if line.startswith('+'):
                # Continuation of previous line
                current_line += ' ' + line[1:].strip()
            else:
                if current_line:
                    processed_lines.append(current_line)
                current_line = line
        if current_line:
            processed_lines.append(current_line)
        
        # Extract MOSFET devices (lines starting with M)
        for line in processed_lines:
            if not line or line.startswith('*') or line.startswith('.'):
                continue
                
            # MOSFET format: Mxxx nd ng ns nb model_name [parameters]
            if line[0].upper() == 'M':
                parts = line.split()
                if len(parts) >= 5:
                    mosfet = {
                        'name': parts[0],
                        'drain': parts[1],
                        'gate': parts[2],
                        'source': parts[3],
                        'bulk': parts[4],
                        'model': parts[5] if len(parts) > 5 else None
                    }
                    mosfets.append(mosfet)
        
        self.mosfets = mosfets
        return mosfets
    
    def get_device_nodes(self, mosfet: Dict) -> Set[str]:
        """Extract non-ground nodes from a MOSFET device."""
        nodes = set()
        
        # Add all nodes except ground (node '0')
        for node_type in ['drain', 'gate', 'source', 'bulk']:
            node = mosfet.get(node_type)
            if node and node != '0':
                nodes.add(node)
        
        return nodes
    
    def build_conflict_graph(self):
        """Build conflict graph where edges connect MOSFETs sharing nodes."""
        n = len(self.mosfets)
        
        # Extract node sets for all MOSFETs
        mosfet_nodes = []
        for mosfet in self.mosfets:
            mosfet_nodes.append(self.get_device_nodes(mosfet))
        
        # Build adjacency list
        self.conflict_graph = defaultdict(set)
        for i in range(n):
            for j in range(i + 1, n):
                # Check if MOSFETs share any nodes
                if mosfet_nodes[i] & mosfet_nodes[j]:  # Set intersection
                    self.conflict_graph[i].add(j)
                    self.conflict_graph[j].add(i)
    
    def greedy_coloring(self) -> int:
        """Apply greedy graph coloring with degree-based ordering."""
        n = len(self.mosfets)
        if n == 0:
            return 0
        
        # Sort vertices by degree (descending) - color high-degree nodes first
        degrees = [(len(self.conflict_graph[i]), i) for i in range(n)]
        degrees.sort(reverse=True)
        
        # Initialize colors
        self.colors = {}
        self.num_colors = 0
        
        # Assign colors
        for _, vertex in degrees:
            # Find the smallest color not used by neighbors
            neighbor_colors = set()
            for neighbor in self.conflict_graph[vertex]:
                if neighbor in self.colors:
                    neighbor_colors.add(self.colors[neighbor])
            
            # Find minimum available color
            color = 0
            while color in neighbor_colors:
                color += 1
            
            self.colors[vertex] = color
            self.num_colors = max(self.num_colors, color + 1)
        
        return self.num_colors
    
    def get_color_groups(self) -> Dict[int, List[int]]:
        """Group MOSFETs by their assigned colors."""
        color_groups = defaultdict(list)
        for vertex, color in self.colors.items():
            color_groups[color].append(vertex)
        return dict(color_groups)
    
    def analyze_netlist(self, netlist_content: str) -> Tuple[int, Dict]:
        """Complete analysis: parse, build graph, color, and return results."""
        # Parse netlist
        self.parse_netlist(netlist_content)
        
        if not self.mosfets:
            return 0, {}
        
        # Build conflict graph
        self.build_conflict_graph()
        
        # Apply coloring
        num_colors = self.greedy_coloring()
        
        # Get color groups
        color_groups = self.get_color_groups()
        
        return num_colors, color_groups
    
    def print_summary(self):
        """Print concise summary of results."""
        print(f"MOSFETs found: {len(self.mosfets)}")
        print(f"Colors needed: {self.num_colors}")
        
        if self.mosfets and self.num_colors > 0:
            # Show color group sizes
            color_groups = self.get_color_groups()
            group_sizes = [len(color_groups[c]) for c in sorted(color_groups.keys())]
            print(f"Group sizes: {group_sizes}")
            
            # Calculate parallelization metrics
            avg_group_size = len(self.mosfets) / self.num_colors
            max_group_size = max(group_sizes)
            efficiency = avg_group_size / max_group_size
            print(f"Parallelization efficiency: {efficiency:.1%}")
    
    def print_detailed(self):
        """Print detailed analysis results."""
        print("\n" + "="*60)
        print("DETAILED ANALYSIS")
        print("="*60)
        
        print(f"\nTotal MOSFETs: {len(self.mosfets)}")
        print(f"Number of colors needed: {self.num_colors}")
        
        if self.mosfets:
            # Calculate graph statistics
            edges = sum(len(neighbors) for neighbors in self.conflict_graph.values()) // 2
            max_degree = max(len(neighbors) for neighbors in self.conflict_graph.values()) if self.conflict_graph else 0
            avg_degree = (2 * edges / len(self.mosfets)) if self.mosfets else 0
            
            print(f"\nConflict Graph Statistics:")
            print(f"  Edges: {edges}")
            print(f"  Maximum degree: {max_degree}")
            print(f"  Average degree: {avg_degree:.2f}")
            
            # Show color groups
            color_groups = self.get_color_groups()
            print(f"\nColor Groups ({self.num_colors} groups):")
            for color in sorted(color_groups.keys()):
                group = color_groups[color]
                mosfet_names = [self.mosfets[i]['name'] for i in group]
                print(f"  Color {color}: {len(group):3d} MOSFETs - {', '.join(mosfet_names)}")
            
            # Show node sharing analysis
            print(f"\nNode Sharing Analysis:")
            node_usage = defaultdict(list)
            for i, mosfet in enumerate(self.mosfets):
                nodes = self.get_device_nodes(mosfet)
                for node in nodes:
                    node_usage[node].append(mosfet['name'])
            
            # Sort nodes by usage count
            sorted_nodes = sorted(node_usage.items(), key=lambda x: len(x[1]), reverse=True)
            print(f"  Most connected nodes (top 10):")
            for node, mosfets in sorted_nodes[:10]:
                print(f"    {node:10s}: {len(mosfets):3d} MOSFETs")


def main():
    """Main function to handle command-line execution."""
    
    # Check command-line arguments
    if len(sys.argv) < 2:
        print("Usage: python3 color.py netlist.cir [options]")
        print("\nOptions:")
        print("  -v, --verbose    Show detailed analysis")
        print("  -q, --quiet      Show only the color count")
        print("\nExample:")
        print("  python3 color.py circuit.cir")
        print("  python3 color.py circuit.cir -v")
        sys.exit(1)
    
    # Parse arguments
    netlist_file = sys.argv[1]
    verbose = False
    quiet = False
    
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            if arg in ['-v', '--verbose']:
                verbose = True
            elif arg in ['-q', '--quiet']:
                quiet = True
    
    # Check if file exists
    if not os.path.exists(netlist_file):
        print(f"Error: File '{netlist_file}' not found.")
        sys.exit(1)
    
    # Read the netlist file
    try:
        with open(netlist_file, 'r') as f:
            netlist_content = f.read()
    except Exception as e:
        print(f"Error reading file '{netlist_file}': {e}")
        sys.exit(1)
    
    # Create analyzer and process netlist
    analyzer = SPICEColoringAnalyzer()
    
    if not quiet:
        print(f"Analyzing: {netlist_file}")
        print("-" * 40)
    
    # Analyze the netlist
    num_colors, color_groups = analyzer.analyze_netlist(netlist_content)
    
    # Output results based on mode
    if quiet:
        # Quiet mode: only output the number
        print(num_colors)
    elif verbose:
        # Verbose mode: detailed analysis
        analyzer.print_summary()
        analyzer.print_detailed()
        print("\n" + "="*60)
        print(f"RESULT: {num_colors} colors needed for parallel execution")
        print("="*60)
    else:
        # Normal mode: summary
        analyzer.print_summary()
        print("-" * 40)
        print(f"Result: {num_colors} colors")
    
    # Return the number of colors as exit code (0 if successful)
    return 0


if __name__ == "__main__":
    sys.exit(main())