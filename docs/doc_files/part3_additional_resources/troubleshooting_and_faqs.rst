.. _troubleshooting:

Troubleshooting and FAQs
============================

Please check out the FAQs below for common issues and solutions. If you encounter a problem that is not listed here, please `let us know. <https://github.com/CCEMGroupTCD/DART/issues>`_

FAQs
-------------------

What's the terminal?
    The terminal is a text-based interface that allows you to interact with your computer using commands. It is also known as the command line or shell. The terminal is a powerful tool that allows you to perform a wide range of tasks, from navigating your file system to executing software like DART.

    How to open your terminal?

    - On Windows, search for Command Prompt or PowerShell in the Start menu.
    - On macOS, open the Applications folder, then Utilities, and double-click on Terminal.
    - On Linux, you can usually find the Terminal in your applications menu or you can use a keyboard shortcut like Ctrl+Alt+T.

Troubleshooting
-------------------

The input yaml file cannot be read:
    Yaml files are quite sensitive to indentation. Make sure that the indentation is correct. You can also use the :ref:`configs module <module_overview>` of DART to get templates for the yaml files.

The input yaml file cannot be found:
    Make sure that the path to the yaml file is correct.

The input yaml file is not valid:
    Verify that the options in the YAML file are set up correctly to match the criteria desired. Ensure that the syntax of the yaml file, including lists and nested elements, is correct. You can use the :ref:`configs module <module_overview>` of DART to get templates for the yaml files.