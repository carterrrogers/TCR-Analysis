# Define the sequences
sequence1 = "IQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDT"
                       
sequence2 = "IQDPDPAVYQLRSPKSSNTSVCLFTDFDSEANVPQSTESTVFSSNSTVLDMRSMDSKSNGALAWSNSTDFKCNSTFNQTFYPTI"

# Compare the two sequences
mismatches = sum(1 for fc, bc in zip(sequence1, sequence2) if fc != bc)
print(f"Total mismatches: {mismatches}")
