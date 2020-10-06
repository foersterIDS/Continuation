Custom documentation
- For a basic insight visit the mathworks page for custom documentation: https://mathworks.com/help/matlab/matlab_prog/display-custom-documentation.html

How to fill the documentation:

To create a documentation i.e. for a function follow these steps:

1. Open MATLAB and create a new live script (.mlx), give it the name of the function or part you'd like to describe
2. Click 'text' and use #title as the title line
3. Use the TEXT and CODE tabs for formatting (see continuation.mlx)
4. Save the file
5. Go to the "Save" button and click the arrow. Use "Export to HTML" to create a html file with the same name as the .mlx file (selected by default)
6. Integrate the html file by opening helptoc.xml using the editor. There are different levels of indent depending on the level of the item in the table of contents.
	First level: main part, here's where Continuation is explained
	Second level: detailed part, here's where functions and example are shown
   Add a new level if necessary by using tab.
   If you want i.e. to include a new documentation for the function 'custom_function' add 
   <tocitem target="custom_function.html">custom_function</tocitem>


How to display the custom documenation:

1. Open MATLAB and go to HOME->Add-Ons->Package Toolbox
2. Click the plus symbol next to "Add toolbox folder" and select the folder of "Continuation". You can also add a version number, an author, an email, a company and give a brief summary and a description
3. Scroll down to "Help Browser Integration" and use "Preview" to see a preview of the documentation.
4. Click Package to create the toolbox
5. Save the file
6. Double click Continuation.mltbx to install the toolbox.
7. Open Help and find your toolbox under "Supplemental Software" at the end of the page.
