
# 🧭 Map Routing Algorithm - C++ Project

## 📌 Overview

This project implements a **shortest-time pathfinding system** across a road network modeled as a graph. It combines both **walking and driving** segments to route from an arbitrary source to an arbitrary destination using an optimized version of **Dijkstra’s Algorithm**.

The goal is to determine the **fastest route** from source to destination, accounting for:

- Variable road lengths and speed limits.
- Free movement (walking) within a radius `R` from source and to destination.
- Real-world map simulation using coordinates.

---

## 📁 Project Structure

```
Map-Routing/
│
├── src/
│   └── main.cpp                 # Main C++ implementation
│
├── output/
│   └── results.txt              # Final output per test case
│
├── README.md                    # This documentation file
│
└── [Test Cases Folder] 🔗       # Downloadable separately (see below)
```

---

## 📦 Features

- Efficient graph construction with adjacency lists.
- Calculates walking paths from source to nearest intersection.
- Calculates driving route using Dijkstra's algorithm.
- Calculates walking segment to destination.
- Precise floating-point output.
- Measures and logs execution time (with/without I/O).

---

## 🧪 Test Cases (Download Required)

This project depends on a **separate folder** of test cases, categorized by size:

| Test Case Type | Max Intersections | Max Roads | Max Queries |
|----------------|-------------------|-----------|-------------|
| Sample         | ≤ 20              | ≤ 50      | ≤ 10        |
| Medium         | ≤ 20,000          | ≤ 25,000  | ≤ 1,000     |
| Large          | ≤ 200,000         | ≤ 250,000 | ≤ 1,000     |

### ⚠️ Required Step:

Before running the program, you **must**:

1. Download the `Test Cases` folder from the provided Google Drive link.
2. Place it in the root directory of the project (or anywhere you prefer).
3. **Update the test case file paths manually** in `main.cpp` between **lines 353 and 373**, as shown below:

```cpp
// Example (inside sample function):
string path_sample = "YOUR/LOCAL/PATH/map.txt";
string path1 = "YOUR/LOCAL/PATH/queries.txt";
```

---

## 🧰 Build & Run Instructions

### 📌 Requirements

- C++17-compatible compiler (e.g., g++, clang++)
- Standard terminal or IDE with file I/O support

### 🔨 Build

```bash
g++ -std=c++17 -O2 -o map_routing src/main.cpp
```

### 🚀 Run

```bash
./map_routing
```

> Ensure your input and output files are placed correctly and file paths are set via `freopen()` or command line.

---

## 🧾 Output Format

For each query, the program prints:

- The list of intersection IDs used in the driving path (excluding virtual walking nodes).
- The total time to reach the destination (in minutes).
- The total path length (in km).
- The total walking distance (in km).
- The total driving distance (in km).
- The total execution time (with/without I/O).

**Example Output:**

```
0 3 4 5 2
4.63 mins
1.72 km
0.28 km
1.44 km
Total I/O Time: 7 ms
Total Execution Time: 13 ms
```

---

## 📚 Implementation Highlights

- **Dijkstra’s Algorithm** is used with a priority queue to ensure optimal time-based routing.
- The graph is constructed dynamically based on the input file.
- Walking zones are modeled by adding virtual edges from/to intersections within the walking radius `R`.

---

## 🧠 Bonuses (Advanced Options)

You can extend the project with:

1. **Speed Intervals**: Variable road speed at different times.
2. **Visualization**: Graphical rendering of map and paths.
3. **Algorithm Optimization**: Sub-logarithmic Dijkstra or A* for better performance.

---

## 📌 Deliverables

- ✅ Fully working C++ codebase
- ✅ Output files for sample/medium/large test cases
- ✅ Code-level documentation and performance analysis

---

## 🏁 Notes

> Remember to:
> - Download the test cases folder.
> - Update file paths from **line 353 to 373** in your code based on the test case folder's URLs.
> - Test each test case by clicking 1 , 2, 3.

---

## 📍 Developed For

**Algorithms Course**  
Department of Information Systems  
Ain Shams University
