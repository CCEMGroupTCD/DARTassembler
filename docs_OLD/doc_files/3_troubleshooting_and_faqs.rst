.. _troubleshooting:

Troubleshooting and FAQs
============================

Please check out this page for common issues and solutions. If you encounter a problem that is not listed here, please `let us know. <https://github.com/CCEMGroupTCD/DART/issues>`_

Issues during installation
----------------------------

Error: libXrender.so.1: cannot open shared object file: No such file or directory
    This error occurs when the required library ``libxrender1`` is missing on your system. On Mac, it is usually preinstalled, but on Linux often not. If you get this issue, you can fix it by installing this library. On Ubuntu for example, you can install it by running the following command:

    .. code-block:: bash

        sudo apt-get install libxrender1

FAQs
-------------------

What's the terminal?
    The terminal is a text-based interface that allows you to interact with your computer using commands. It is also known as the command line or shell. The terminal is a powerful tool that allows you to perform a wide range of tasks, from navigating your file system to executing software like DART.

    How to open your terminal?

    - On Windows, search for Command Prompt or PowerShell in the Start menu.
    - On macOS, open the Applications folder, then Utilities, and double-click on Terminal.
    - On Linux, you can usually find the Terminal in your applications menu or you can use a keyboard shortcut like Ctrl+Alt+T.

What is pip?
    pip is a widely used package manager for Python that allows you to install and manage Python packages. You can use pip to install DART and its dependencies.

What is a yaml file?
    YAML is a human-readable file format that is commonly used for configuration files. YAML files are easy to read and write, making them a popular choice for storing settings and options for software applications.

How to edit a yaml file?
    You can edit a yaml file using a text editor like Notepad, TextEdit, or Visual Studio Code. Simply open the yaml file in your text editor, make the desired changes, and save the file. Be careful with the indentation in yaml files, as it is significant and can affect the structure of the file.

Troubleshooting
-------------------

The input yaml file cannot be found:
    Make sure that the path to the yaml file is correct.

The input yaml file is not valid or cannot be read:
    Yaml files are quite sensitive to indentation. Make sure that the indentation is correct. Verify that the options in the YAML file are set up correctly to match the criteria desired. In case of any issues with the yaml configuration files, you can use the :ref:`configs module <module_overview>` of DART to get two templates files so that you can fill in the required information without worrying about the syntax. You can also use tools like the `online yaml linter <https://www.yamllint.com/>`_ to make sure that your yaml file is valid.

