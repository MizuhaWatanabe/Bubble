#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define max_data 10000

int main(){

    FILE *input;
    char row_data[max_data];
    char *name_key, *name_value,
         *univ_key, *univ_value,
         *tel_name, *tel_value,
         *height_name, *height_value,
         *weight_name, *weight_value;
    char *dummy_index;
    char skip[] = ":,";

    input = fopen("test.csv", "r");
    fscanf(input, "%s", row_data);

    name_key = strtok(row_data, skip);
    printf("%s\n", name_key);

    while(1){
        name_value = strtok(NULL, skip);
        dummy_index = strtok(NULL, skip);
        if(dummy_index != NULL){
            name_key = dummy_index;
        }
        else{
            break;
        }

        printf("%s,%s\n", name_key, name_value);
    }



}
