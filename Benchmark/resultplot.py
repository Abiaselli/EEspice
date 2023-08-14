import pandas as pd
import matplotlib.pyplot as plt

def plot_data(df, algorithm_name):
    df = df[df["Factorization"] != -1]
    
    thread_counts = df["threads"].unique()
    for thread in thread_counts:
        filtered_data = df[df["threads"] == thread]
        
        # Sort the filtered data by the 'nnz' column
        filtered_data = filtered_data.sort_values(by="rows")
        plt.loglog(filtered_data["rows"], filtered_data["Factorization"], label=f'Algorithm={algorithm_name}, Threads={thread}')
        # plt.plot(filtered_data["rows"], filtered_data["Factorization"], label=f'Algorithm={algorithm_name}, Threads={thread}')

# Load data from both CSV files
df1 = pd.read_csv("results_nicslu_kernel.csv")
# df2 = pd.read_csv("results_klu_kernel.csv")

# Assuming that the 'algorithmname' column exists and contains the name of the algorithm
plot_data(df1, df1['algorithmname'].iloc[0])
# plot_data(df2, df2['algorithmname'].iloc[0])

plt.xlabel('Matrix Size')
plt.ylabel('Time for Factorization')
plt.title('Factorization Time vs Number of non-zeros')
plt.legend()
plt.grid(True)
plt.tight_layout()

# plt.show()
plt.savefig('resultplot.png',dpi = 600)


