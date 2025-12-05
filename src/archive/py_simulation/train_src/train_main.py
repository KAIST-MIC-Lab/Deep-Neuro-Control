import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import torch
import math
import pickle
import numpy as np
import matplotlib.pyplot as plt
from SimpleNN import SimpleNN

# DATA_NAME = "src/py_simulation/train_src/result.pickle"
DATA_NAME = "result.pickle"

EPOCH = 100
REPORT_EPOCH = 1
LEARNING_RATE = 1e-1
BATCH_SIZE = 128
MINI_BATCH_SIZE = 32

VAL_RATIO = 0.2

def data_load():
    data = pickle.load(open(DATA_NAME, "rb"))
    keys = list(data.keys())

    x_hist =    np.array([])
    r_hist =    np.array([])
    ud_hist =   np.array([])
    M_hist =    np.array([])

    for key_idx in range(len(keys)):
        key = keys[key_idx]
    
        x_hist =    np.append(x_hist,   (data[key]["x_hist"] ))
        r_hist =    np.append(r_hist,   (data[key]["r_hist"] ))
        ud_hist =   np.append(ud_hist,  (data[key]["ud_hist"]))
        M_hist =    np.append(M_hist,   (data[key]["M_hist"] ))

    x_hist =    x_hist.reshape(-1, 2,1)
    r_hist =    r_hist.reshape(-1, 2,1)
    ud_hist =   ud_hist.reshape(-1, 2,1)
    M_hist =    np.linalg.cholesky(M_hist.reshape(-1, 2,2))

    M_hist = M_hist.reshape(-1, 4, 1)
    M_hist = M_hist[:,[0,2,3], :]

    # nn_inputs = np.concatenate((x_hist, r_hist, ud_hist), axis=1)
    nn_inputs = np.concatenate((x_hist, r_hist), axis=1)
    nn_outputs = M_hist

    # Shuffle the data
    indices = np.arange(len(nn_inputs))
    np.random.shuffle(indices)
    nn_inputs = nn_inputs[indices]
    nn_outputs = nn_outputs[indices]

    # Split the data into training and validation sets
    split = int(len(nn_inputs) * (1 - VAL_RATIO))
    train_inputs = nn_inputs[:split]
    train_outputs = nn_outputs[:split]
    val_inputs = nn_inputs[split:]
    val_outputs = nn_outputs[split:]

    print(f"Train inputs shape: {train_inputs.shape}")
    print(f"Train outputs shape: {train_outputs.shape}")

    return train_inputs, train_outputs, val_inputs, val_outputs

def train(train_inputs, train_outputs, val_inputs, val_outputs):
    train_result = []

    dtype = torch.float
    device = torch.device("mps") if torch.backends.mps.is_available() else torch.device("cpu")
    print(f"Using {device} device")

    nn = SimpleNN()    
    nn = nn.to(device)
    loss_fn = torch.nn.MSELoss(reduction='mean')

    BATCH_NUM = int(len(train_inputs) / BATCH_SIZE)
    MINI_BATCH_NUM = int(BATCH_SIZE / MINI_BATCH_SIZE)

    for epoch_idx in range(EPOCH):
        print(f"Epoch {epoch_idx + 1}/{EPOCH}")

        for batch_idx in range(BATCH_NUM):

            batch_inputs = train_inputs[batch_idx:batch_idx + BATCH_SIZE]
            batch_outputs = train_outputs[batch_idx:batch_idx + BATCH_SIZE]

            for mini_batch_idx in range(MINI_BATCH_NUM):
                print(f"Batch {batch_idx + 1}/{BATCH_NUM} Mini Batch {mini_batch_idx + 1}/{MINI_BATCH_NUM}", end="\r")

                mini_batch_inputs = batch_inputs[mini_batch_idx:mini_batch_idx + MINI_BATCH_SIZE]
                mini_batch_outputs = batch_outputs[mini_batch_idx:mini_batch_idx + MINI_BATCH_SIZE]

                x = torch.tensor(mini_batch_inputs, device=device, dtype=dtype)
                y = torch.tensor(mini_batch_outputs, device=device, dtype=dtype)

                x = x.reshape(-1, 4)
                y = y.reshape(-1, 3)

                # optimizer = torch.optim.RMSprop(nn.parameters(), lr=LEARNING_RATE)
                optimizer = torch.optim.SGD(nn.parameters(), lr=LEARNING_RATE)
                for epoch_idx in range(EPOCH):
                    # Forward pass: compute predicted y by passing x to the nn.
                    y_pred = nn(x)

                    # Compute and print loss.
                    loss = loss_fn(y_pred, y)

                    optimizer.zero_grad()
                    loss.backward()
                    optimizer.step()

        # Validation
        with torch.no_grad():
            val_idx = np.random.randint(0, len(val_inputs), BATCH_SIZE)
            
            val_x = torch.tensor(val_inputs[val_idx], device=device, dtype=dtype)
            val_y = torch.tensor(val_outputs[val_idx], device=device, dtype=dtype)

            val_x = val_x.reshape(-1, 4)
            val_y = val_y.reshape(-1, 3)

            val_y_pred = nn(val_x)
            val_loss = loss_fn(val_y_pred, val_y)
            print(f"Validation Loss: {val_loss.item()}")
            train_result.append(val_loss.item())

    return nn, train_result

if __name__ == "__main__":

    train_inputs, train_outputs, val_inputs, val_outputs = data_load()
    print("Data loaded")

    nn, train_result = train(train_inputs, train_outputs, val_inputs, val_outputs)

    torch.save(nn.state_dict(), "model.pth")
    print("Model saved")
    print("Training result:")
    print(train_result)
    plt.plot(train_result)
