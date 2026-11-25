# Mirror Furnace Ray Tracer

Explore how light travels through a two-mirror furnace. Move the source and mirrors, adjust their tilt, and see where rays cross a horizontal detector. The simulator runs on your computer and opens in your web browser. It uses Python's standard library, so there are no extra Python packages to install.

## 1. Install Python

Use Python **3.9 or later**. If you're installing it for the first time, choose the latest stable Python 3 release for your computer.

### Windows

1. Open the [Python for Windows downloads page](https://www.python.org/downloads/windows/) and install the **Python install manager**. If you use the traditional installer instead, check **Add python.exe to PATH** when the installer offers it.
2. Open the Start menu, type **PowerShell**, and open it. PowerShell is where you type the commands below. Press **Enter** after each command.
3. Check that Python is available.

   ```powershell
   python --version
   ```

   The install manager may install its default Python version the first time you run this command. You should then see a version such as `Python 3.14.6`. If PowerShell can't find `python` but `py --version` works, use `py` in place of `python` in the Windows commands below. If neither works, close and reopen PowerShell after installing Python.

### macOS

1. Open the [Python for macOS downloads page](https://www.python.org/downloads/macos/) and download the **macOS installer** for the latest stable Python 3 release. Double-click the downloaded `.pkg` file and follow the installer.
2. Press **Command–Space**, type **Terminal**, and open it. Type the command below and press **Return**.

   ```bash
   python3 --version
   ```

   You should see a version such as `Python 3.14.6`.

### Linux

1. Open **Terminal** from your applications menu and check whether Python is already installed.

   ```bash
   python3 --version
   ```

2. If you see Python 3.9 or later, continue to the next section. Otherwise, install Python 3 with your distribution's package manager.

   **Ubuntu or Debian**

   ```bash
   sudo apt update
   sudo apt install python3
   ```

   **Fedora**

   ```bash
   sudo dnf install python3
   ```

   Run `python3 --version` again to confirm. If your distribution provides a version older than 3.9, use its instructions to install a newer Python 3.

## 2. Download and extract the simulator

1. At the top of this GitHub repository, click **Code**, then **Download ZIP**. You can also [download the ZIP directly](https://github.com/Plexanovmax/mirror-furnace-ray-tracer/archive/refs/heads/main.zip).
2. Open your **Downloads** folder and extract the ZIP. On Windows, right-click the ZIP and choose **Extract All**. On macOS, double-click the ZIP. On Linux, right-click the ZIP and choose **Extract Here**, or use your archive manager.
3. Open the extracted folder that contains `mirror_furnace.py` and `furnace_drawing.png`. Keep both files in that folder so the simulator can load the drawing.

## 3. Start the simulator

### Windows

In File Explorer, open the extracted folder containing `mirror_furnace.py`. Click the address bar, type `powershell`, and press **Enter** to open PowerShell in that folder. Then run the following command.

```powershell
python .\mirror_furnace.py
```

### macOS or Linux

Open Terminal. If you extracted the ZIP in **Downloads**, run these commands one line at a time.

```bash
cd "$HOME/Downloads/mirror-furnace-ray-tracer-main"
python3 mirror_furnace.py
```

If you extracted the ZIP somewhere else, change the path after `cd` to match that folder. On macOS, you can type `cd ` with a space after it, drag the extracted folder from Finder into Terminal to insert its path, and press **Return**.

### Open or stop the simulator

The simulator normally opens your default browser automatically. If it doesn't, copy the address printed in PowerShell or Terminal into your browser's address bar. It starts with `http://127.0.0.1:` and runs on your own computer. Leave the command window open while you use the simulator. To stop it, return to that window and press **Ctrl+C**.

If you see **can't open file**, use `dir` in PowerShell or `ls` in Terminal to check that `mirror_furnace.py` is in the current folder.

## Use the simulator

- Choose **Bottom**, **Top**, **Both**, or **None** under **Point source**. The slider controls the total number of rays, from 2 to 100. **Both** splits that total between the two sources.
- Adjust the source and mirrors using the sliders, number fields, or **+ / −** buttons. **Fine steps** changes the buttons from 1 mm to 0.1 mm. Mirror tilt is specified as rise on a 100 mm reference lever.
- Enable **XY detector plane · 3D** to trace rays in three dimensions. The main drawing shows the XZ projection, and the panel on the right shows crossings on the detector disk. The disk radius is adjustable from 30 to 550 mm, and the rings are 10 mm apart.
- Hover over the furnace to inspect a region in the zoom panel. A smaller lens gives more magnification. Use **PDF opacity** to adjust the background drawing.

The background PNG is a transparent crop of PSI furnace drawing **DRW-00322318**. You can replace `furnace_drawing.png` with an annotated copy using the same crop and alignment, then refresh the browser. The optical model uses ideal truncated mirror surfaces and omits shutters and other mechanical obstructions.
