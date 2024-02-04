import pandas as pd
from adjacency_list_graph import AdjacencyListGraph
from dijkstra import dijkstra
from bellman_ford import bellman_ford
import matplotlib.pyplot as plt
from itertools import combinations


# Creating the Graph
def build_graph(data):

    # Creating a Station to Index mapping dictionary for calling reference
    station_to_index = {station: index for index, station in enumerate(set(data['start']) | set(data['end']))}

    combined_graph = AdjacencyListGraph(card_V=len(station_to_index), directed=True, weighted=True)
    inserted_edges = set()

    # Iterating through the Excel File to add data to the Graph
    for _, row in data.iterrows():
        start_index = station_to_index[row['start']]
        end_index = station_to_index[row['end']]
        duration = row['duration']

        # Insert the original edge
        if (start_index, end_index) not in inserted_edges:
            combined_graph.insert_edge(start_index, end_index, duration)
            inserted_edges.add((start_index, end_index))

        # Insert the reversed edge
        if (end_index, start_index) not in inserted_edges:
            combined_graph.insert_edge(end_index, start_index, duration)
            inserted_edges.add((end_index, start_index))

    return combined_graph, station_to_index

# Djikstra algorithm
def find_shortest_path(graph, start, end):
    _, pi = dijkstra(graph, start)
    shortest_path = [end]
    current = end
    while pi[current] is not None:
        shortest_path.insert(0, pi[current])
        current = pi[current]
    return shortest_path


# Task 1A - Calculating total duration in minutes of the journey path
def calculate_total_duration(graph, path):
    if path:
        # Calculates the total duration based on the durations of edges in the path.
        total_duration = sum(graph.find_edge(path[i], path[i + 1]).get_weight() for i in range(len(path) - 1))
        return total_duration
    else:
        return None


# Task 1A - Using Dijkstra Algorithm to find the shortest journey
def run_task_1(graph, station_to_index):
    print("                    London Underground Route System")
    print("--------------------------------------------------------------------------")
    print("Task 1(a) - using Dijkstra's algorithm.")

    # Taking the input from the user
    start_station = input("Enter the starting station: ")
    end_station = input("Enter the destination station: ")


    # Getting the index value from the mapping dictionary
    start_index = station_to_index.get(start_station)
    end_index = station_to_index.get(end_station)

    # Calculating the shortest path and total duration
    if start_index is not None and end_index is not None:
        shortest_path = find_shortest_path(graph, start_index, end_index)
        total_duration = calculate_total_duration(graph, shortest_path)
        route = [list(station_to_index.keys())[i] for i in shortest_path]
        print("\nRoute Information:")
        if shortest_path:
            print(f"Shortest path: {' -> '.join(route)}")
            print(f"Total Duration (in minutes): {total_duration}")
        else:
            print(f"No path found from {start_station} to {end_station}")
    else:
        print("Invalid station names.")

    return shortest_path, start_station, end_station, total_duration


# TASK 1B - Calculating all the journey times for the whole network to plot the histogram
def calculate_all_journey_times(graph, station_to_index):

    # Creating an empty list to store the journey times
    all_journey_times = []

    # Getting the station name list from the mapping dictionary
    stations = list(station_to_index.keys())

    # Looping until we use the Dijkstra Algorithm on all instances of the stations
    for i in range(len(stations)):
        for j in range(i + 1, len(stations)):
            start_index = station_to_index[stations[i]]
            end_index = station_to_index[stations[j]]

            # Using the Dijkstra Algorithm to find the Journey times
            _, pi = dijkstra(graph, start_index)
            shortest_path = find_shortest_path(graph, start_index, end_index)
            total_duration = calculate_total_duration(graph, shortest_path)

            # Appending the journey times to the list we created
            all_journey_times.append(total_duration)

    # Returns a list of all journey times in the whole network.
    return all_journey_times


# Task 1B - Plotting the Histogram for Journey times between Station Pairs in terms of minutes
def plot_histogram(journey_times):
    plt.hist(journey_times, bins=20, edgecolor='black')
    plt.title('Journey Time Distribution between Station Pairs')
    plt.xlabel('Journey Time (minutes)')
    plt.ylabel('Frequency')
    plt.show()


# Task 2A - Calculating the total stops on the journey path
def calculate_total_stops(path):
    return len(path) - 1 if path else None


# Task 2A - Using the Djikstra Algorithm to find the number of stops between two given stations.
def run_task_2(graph, station_to_index):
    print("\nTask 2(a) - using Dijkstra's algorithm.")

    # Taking input from user
    start_station_task2 = input("Enter the starting station: ")
    end_station_task2 = input("Enter the destination station: ")

    # Getting the index value from the mapping dictionary
    start_index_task2 = station_to_index[start_station_task2]
    end_index_task2 = station_to_index[end_station_task2]

    # Calculating the shortest path and total duration
    shortest_path_task2 = find_shortest_path(graph, start_index_task2, end_index_task2)
    total_duration_task2 = calculate_total_duration(graph, shortest_path_task2)
    index_to_station = {index: station for station, index in station_to_index.items()}

    print("\nRoute Information:")
    if shortest_path_task2:
        station_path_task2 = [index_to_station[i] for i in shortest_path_task2]
        print(f"Route taken: {' -> '.join(station_path_task2)}")
        print(f"Count of Stations or Stops: {calculate_total_stops(shortest_path_task2)}")
    else:
        print(f"No path found from {start_station_task2} to {end_station_task2}")


# TASK 2B - Calculating all the journey stops for the whole network to plot the histogram
def calculate_all_num_stops(graph, station_to_index):
    # Creating an empty list to store the journey times
    all_num_stops = []

    # Getting the station name list from the mapping dictionary
    stations = list(station_to_index.keys())

    # Looping until we use the Dijkstra Algorithm on all instances of the stations
    for i in range(len(stations)):
        for j in range(i + 1, len(stations)):
            start_index = station_to_index[stations[i]]
            end_index = station_to_index[stations[j]]

            # Using the Dijkstra Algorithm to find the Journey stops
            _, pi = dijkstra(graph, start_index)
            shortest_path = find_shortest_path(graph, start_index, end_index)
            num_stops = calculate_total_stops(shortest_path)

            # Appending the journey times to the list we created
            all_num_stops.append(num_stops)

    # Returns a list of all journey stops in the whole network.
    return all_num_stops


# Task 2B - Plotting the Histogram for Journey times between Station Pairs in terms of stops
def plot_histogram_num_stops(num_stops):
    plt.hist(num_stops, bins=20, edgecolor='black')
    plt.title('Number of Stops Distribution between Station Pairs using Bellman Ford Algorithm')
    plt.xlabel('Number of Stops')
    plt.ylabel('Frequency')
    plt.show()


# Task 3A - Using the Bellman Ford Algorithm to find the number of stops between two station pairs
def run_task_3():
    print("\nTask 3(a) - using Bellman Ford algorithm.")

    # Taking input from user
    start_station_task3 = input("Enter the starting station : ")
    end_station_task3 = input("Enter the destination station : ")

    # Map station names to integer indices
    start_index_task3 = station_to_index[start_station_task3]
    end_index_task3 = station_to_index[end_station_task3]

    # Use Bellman-Ford algorithm
    distances, predecessors, no_negative_cycle = bellman_ford(combined_graph, start_index_task3)

    if no_negative_cycle:
        # Extract the shortest path
        shortest_path_task3 = [end_index_task3]
        current = end_index_task3
        while predecessors[current] is not None:
            shortest_path_task3.insert(0, predecessors[current])
            current = predecessors[current]

        # Calculating the Total Duration in Task 3
        total_duration_task3 = calculate_total_duration(combined_graph, shortest_path_task3)

        print("\nRoute Information:")
        route2 = [list(station_to_index.keys())[i] for i in shortest_path_task3]
        print(f"Route taken: {' -> '.join(route2)}")
        print(f"Count of Stations or Stops: {calculate_total_stops(shortest_path_task3)}")
        print(f"Total Duration (in minutes): {total_duration_task3}")
    else:
        print("Negative-weight cycle detected. Cannot find the shortest path.")


# TASK 3B - Calculating all the journey stops using Bellman Ford ALgorithm for the whole network to plot the histogram
def calculate_all_num_stops_bf(graph, station_to_index):

    # Creating an empty list to store the journey times
    all_num_stops_bf = []

    # Getting the station name list from the mapping dictionary
    stations = list(station_to_index.keys())

    # Looping until we use the Bellman Ford Algorithm on all instances of the stations
    for i in range(len(stations)):
        for j in range(i + 1, len(stations)):
            start_index = station_to_index[stations[i]]
            end_index = station_to_index[stations[j]]

            # Using the Bellman Ford Algorithm to find the Journey stops
            _, predecessors, no_negative_cycle = bellman_ford(graph, start_index)

            # Appending the journey stops to the list we created
            if no_negative_cycle:
                shortest_path = [end_index]
                current = end_index
                while predecessors[current] is not None:
                    shortest_path.insert(0, predecessors[current])
                    current = predecessors[current]

                num_stops = calculate_total_stops(shortest_path)
                all_num_stops_bf.append(num_stops)
            else:
                all_num_stops_bf.append(None)

    # Returns a list of all journey stops in the whole network.
    return all_num_stops_bf

# Task 3B - Plotting the Histogram for Journey times between Station Pairs in terms of stops
def plot_histogram_num_stops_bf(num_stops_bf):
    plt.hist(num_stops_bf, bins=20, edgecolor='black')
    plt.title('Number of Stops Distribution between Station Pairs (Bellman-Ford)')
    plt.xlabel('Number of Stops')
    plt.ylabel('Frequency')
    plt.show()


# Task 4A - Finding stations that we can close
def find_viable_adjacent_station_pairs(graph, station_to_index):
    viable_adjacent_pairs = set()

    for _, row in data.iterrows():
        start_index = station_to_index[row['start']]
        end_index = station_to_index[row['end']]
        duration = row['duration']

        # Temporarily remove the edge and check if there is still a path between the stations
        graph.delete_edge(start_index, end_index)
        _, pi = dijkstra(graph, start_index)

        if pi[end_index] is not None:
            # Add the station pair to the set
            pair = tuple(sorted([row['start'], row['end']]))
            viable_adjacent_pairs.add(pair)

        # Add the edge back for the next iteration
        graph.insert_edge(start_index, end_index, duration)

    return viable_adjacent_pairs


# Task 4B - Calculating both journey times in terms of mins and stops
def calculate_all_journey_info(graph, station_to_index):
    # Creating an empty list to store the journey times
    all_journey_info = []

    # Looping until we use the Dijkstra Algorithm on all instances of the stations
    for i in range(len(stations)):
        for j in range(i + 1, len(stations)):
            start_index = station_to_index[stations[i]]
            end_index = station_to_index[stations[j]]

            # Using the Dijkstra Algorithm to find the Journey stops
            _, pi = dijkstra(graph, start_index)
            shortest_path = find_shortest_path(graph, start_index, end_index)
            total_duration = calculate_total_duration(graph, shortest_path)
            num_stops = calculate_total_stops(shortest_path)

            # Appending the journey mins and stops to the list we created
            all_journey_info.append((total_duration, num_stops))

    # Returns a list of all journey mins and stops in the whole network.
    return all_journey_info


# Task 4B - Plotting the histogram for Journey times between Station Pairs in terms of stops and mins
def plot_histogram_journey_info(journey_info):
    total_durations, num_stops = zip(*journey_info)

    bar_width = 0.35
    index = range(len(journey_info))

    fig, ax = plt.subplots()
    bars1 = ax.bar(index, total_durations, bar_width, label='Journey Time (mins)')
    bars2 = ax.bar([i + bar_width for i in index], num_stops, bar_width, label='Number of Stops')

    ax.set_xlabel('Station Pairs')
    ax.set_ylabel('Values')
    ax.set_title('Journey Time and Number of Stops Distribution between Station Pairs')
    ax.set_xticks([i + bar_width / 2 for i in index])
    ax.set_xticklabels([f'{stations[i]} - {stations[j]}' for i, j in combinations(range(len(stations)), 2)], rotation=45, ha='right')
    ax.legend()

    plt.show()


# Reading the Data from the Excel File
def read_data(file_path):
    ds = pd.read_excel(file_path, header=None, names=['line', 'start', 'end', 'duration'])
    ds = ds.dropna() # Deleting the NAN values from the file
    return ds


# File path to your Excel data
file_path = "London Underground data.xlsx"

# Read data from Excel
data = read_data(file_path)

# Build the graph
combined_graph, station_to_index = build_graph(data)

# Getting the station list from the mapping dictionary
stations = list(station_to_index.keys())


# Run the tasks

# TASK 1A
run_task_1(combined_graph, station_to_index)

print("--------------------------------------------------------------------------")

print("Generating Histogram for TASK 1B")

print("--------------------------------------------------------------------------")


# TASK 2A
run_task_2(combined_graph, station_to_index)

print("--------------------------------------------------------------------------")

print("Generating Histogram for TASK 2B")

print("--------------------------------------------------------------------------")

# TASK 3A
run_task_3()

print("--------------------------------------------------------------------------")

print("Generating Histogram for TASK 3B")

print("--------------------------------------------------------------------------")


# Task 4A
# Find the list of adjacent station pairs where travel remains viable without modifying the original graph
viable_adjacent_station_pairs = find_viable_adjacent_station_pairs(combined_graph, station_to_index)

# Print the unique list of adjacent station pairs where travel remains viable even if the connection is severed
print("Task 4(a)- List of adjacent station pairs where travel remains viable even if the connection is severed:")
for pair in viable_adjacent_station_pairs:
    print(pair[0], "-", pair[1])

# TASK 2B - Calculate all number of stops
all_num_stops = calculate_all_num_stops(combined_graph, station_to_index)

# TASK 2B - Plot the histogram for number of stops
plot_histogram_num_stops(all_num_stops)


# TASK 3B - Calculate all number of stops using Bellman-Ford
all_num_stops_bf = calculate_all_num_stops_bf(combined_graph, station_to_index)

# TASK 3B - Plot the histogram for number of stops using Bellman-Ford
plot_histogram_num_stops_bf(all_num_stops_bf)

# TASK 4B - Calculate all journey information
all_journey_info = calculate_all_journey_info(combined_graph, station_to_index)

# TASK 4B - Plot the grouped bar chart for journey time and number of stops
plot_histogram_journey_info(all_journey_info)