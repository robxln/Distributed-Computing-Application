# Distributed Computing Application

## Overview

This README provides a comprehensive guide to a distributed computing application developed to handle large data processing tasks using varying communication strategies based on the `commError` parameter. This document will give you an understanding of how the code is structured, how it functions, and how to run it effectively.

## Table of Contents

- [Overview](#overview)
- [Implementation Details](#implementation-details)
  - [Communication Strategies](#communication-strategies)
  - [Topology Discovery](#topology-discovery)
  - [Calculations](#calculations)

## Overview

The distributed computing application aims to process a large data array efficiently by distributing the workload across a cluster of workers. The key feature of the application is its ability to adapt to two primary communication strategies based on the `commError` parameter:

- **commError == 0:** Utilizes a ring communication algorithm for data transfer and collection.
- **commError == 1/2:** Implements a communication strategy similar to Depth-First Search (DFS) for data transfer.

## Implementation Details

### Communication Strategies

#### Ring Communication (commError == 0)

In this mode, the application employs a ring communication strategy for data transmission and collection. The process unfolds as follows:

1. All clusters must synchronize and discover the complete system topology.
2. The topology information is transmitted between clusters through a "ring" structure.
3. Once it is ensured that all clusters have the correct topology information, the ring communication begins, passing data sequentially from one cluster to another until it returns to the source cluster.

#### DFS-Like Communication (commError == 1/2)

For `commError` values of 1 or 2, the application follows a communication strategy inspired by Depth-First Search. The following steps outline this approach:

1. Each cluster initiates a DFS-like traversal to determine the topology and communication paths.
2. Topology information is transmitted incrementally, with each cluster forwarding it to the next one along the path.
3. The communication path spans all clusters in the system to ensure efficient data distribution to every worker.

### Topology Discovery

In both communication strategies, discovering the system's topology is crucial. It is achieved through a parent vector where, for each cluster, the parent cluster is identified. This vector aids in organizing communication paths effectively. The synchronization and topology discovery process may vary depending on the communication strategy used.

### Calculations

The primary objective of the application is to process a large data array in parallel. The workload distribution works as follows:

- Process 0 allocates workload to available workers, ensuring balanced work distribution based on the array's size and the number of workers.
- Workload allocation is optimized for minimal processing time.

The computation process follows these steps:

1. The source cluster, typically Cluster 0, allocates the workload and initiates data transfer.
2. Each cluster processes its assigned portion of the array (multiplying elements by 5).
3. Partial results are communicated sequentially along the designated communication path.
4. The final result is accumulated and displayed.
