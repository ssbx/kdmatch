/*
 * Load and manipulate and perform actions on Sample common to a
 * single catalog.
 *
 * Copyright (C) 2017 University of Bordeaux. All right reserved.
 * Written by Emmanuel Bertin
 * Written by Sebastien Serre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef __CATALOG_H__
#define __CATALOG_H__

#include "scamp.h"

/**
 * Open a catalog. Presently only support sextractor catalogs. The Field
 * structure given in input must be freed by the user with Catalog_free().
 *
 * Tread safe.
 */
extern void
Catalog_open(char *file, Field *field);

/**
 * Print the content of catalogs. Used for debugging purpose.
 *
 * Thread safe.
 */
extern void
Catalog_dump(Field *field);

/**
 * Free all memory allocated for a field.
 *
 * Thread safe.
 */
extern void
Catalog_freeField(Field *field);



#endif /* __CATALOG_H__ */
